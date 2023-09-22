process FILTER_AND_FIX_STRAND_FLIPS {

    input:
    tuple val(filename), path(map_file)
    path(strand_file)
    path(refalt_file)

    output:
    path "${filename}.vcf.gz", emit: vcf_files
    path "${filename}.vcf.gz.tbi", emit: vcf_files_index
    path "${filename}.qc.samples", emit: samples_runs
    path "${filename}.qc.snps", emit: snps_runs
    path "*.statistics", emit: filter_statistics

    script:
    """

    # count all lines in map file and ped file as write as step0 to .statistics
    total_samples=\$(cat ${filename}.ped | wc -l)
    total_snps=\$(cat ${filename}.map | wc -l)
    printf "name samples snps\n" > ${filename}.statistics
    printf "step00 \${total_samples} \${total_snps}" >> ${filename}.statistics

    # Step 1: Remove all indels, "I" and "D"
    plink --file ${filename} \
        --snps-only just-acgt \
        --make-bed \
        --out ${filename}.step01

    plink-statistics "step01" ${filename}.step01 ${filename}.statistics


    # Step 2: Remove all non-autosomale SNPs
    plink --bfile ${filename}.step01 \
        --chr 1-22 \
        --make-bed \
        --out ${filename}.step02

    plink-statistics "step02" ${filename}.step02 ${filename}.statistics


    # Step 3: Update strand flips (https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh)
    update_build.sh \
        ${filename}.step02 \
        ${strand_file} \
        ${filename}.step03

    plink-statistics "step03" ${filename}.step03 ${filename}.statistics

    # Step 4: Remove all autosomale snps after update strand flips
    plink --bfile ${filename}.step03 \
        --chr 1-22 \
        --make-bed \
        --out ${filename}.step04

    plink-statistics "step04" ${filename}.step04 ${filename}.statistics


    # Step 5: Harmonize ref/alt alleles and retain only SNPs in the refalt file.
    plink --bfile ${filename}.step04 \
        --extract ${refalt_file} \
        --a2-allele ${refalt_file} \
        --recode ${params.useDoubleId ? 'vcf-iid' : 'vcf'} \
        --out ${filename}.harmonized

    # remove all variants that have no ref allele inside

    grep  "Warning: Impossible A2 allele assignment for variant *" ${filename}.harmonized.log | awk '{print substr(\$NF,1,length(\$NF)-1)}' > ${filename}.harmonized.snps

    vcftools --vcf ${filename}.harmonized.vcf \
        --exclude ${filename}.harmonized.snps \
        --recode --stdout | bgzip -c > ${filename}.vcf.gz

    tabix ${filename}.vcf.gz

    java -jar /opt/genomic-utils.jar vcf-statistics --name "step05" --input "${filename}.vcf.gz" --output "${filename}.statistics"

    # Calculate snp call rate and sample call rate (per run)
    java -jar /opt/genomic-utils.jar vcf-quality-control \
        ${filename}.vcf.gz \
        --minSnpCallRate ${params.minSnpCallRate}  \
        --minSampleCallRate ${params.minSampleCallRate}  \
        --chunkSize ${params.chunkSize} \
        --output ${filename}.qc
    """

}