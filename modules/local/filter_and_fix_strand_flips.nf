process FILTER_AND_FIX_STRAND_FLIPS {

    input:
    tuple val(run), val(filename), path(map_file), path(ped_file)
    path(strand_file)
    path(refalt_file)

    output:
    path "${run}.vcf.gz", emit: vcf_files
    path "${run}.vcf.gz.tbi", emit: vcf_files_index
    path "${run}.qc.samples", emit: samples_runs
    path "${run}.qc.snps", emit: snps_runs
    path "*.statistics", emit: filter_statistics

    script:
    """

    # count all lines in map file and ped file as write as step0 to .statistics
    total_samples=\$(cat ${filename}.ped | wc -l)
    total_snps=\$(cat ${filename}.map | wc -l)
    printf "name samples snps\n" > ${run}.statistics
    printf "step00 \${total_samples} \${total_snps}" >> ${run}.statistics

    # Step 1: Remove all indels, "I" and "D"
    plink --file ${filename} \
        --snps-only just-acgt \
        --make-bed \
        --out ${run}.step01

    plink-statistics "step01" ${run}.step01 ${run}.statistics


    # Step 2: Remove all non-autosomale SNPs
    plink --bfile ${run}.step01 \
        --chr 1-22 \
        --make-bed \
        --out ${run}.step02

    plink-statistics "step02" ${run}.step02 ${run}.statistics


    # Step 3: Update strand flips (https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh)
    update_build.sh \
        ${run}.step02 \
        ${strand_file} \
        ${run}.step03

    plink-statistics "step03" ${run}.step03 ${run}.statistics

    # Step 4: Remove all autosomale snps after update strand flips
    plink --bfile ${run}.step03 \
        --chr 1-22 \
        --make-bed \
        --out ${run}.step04

    plink-statistics "step04" ${run}.step04 ${run}.statistics


    # Step 5: Harmonize ref/alt alleles and retain only SNPs in the refalt file.
    plink --bfile ${run}.step04 \
        --extract ${refalt_file} \
        --a2-allele ${refalt_file} \
        --recode ${params.useDoubleId ? 'vcf-iid' : 'vcf'} \
        --out ${run}.harmonized

    # remove all variants that have no ref allele inside

    grep  "Warning: Impossible A2 allele assignment for variant *" ${run}.harmonized.log | awk '{print substr(\$NF,1,length(\$NF)-1)}' > ${run}.harmonized.snps

    vcftools --vcf ${run}.harmonized.vcf \
        --exclude ${run}.harmonized.snps \
        --recode --stdout | bgzip -c > ${run}.vcf.gz

    tabix ${run}.vcf.gz

    java -jar /opt/genomic-utils.jar vcf-statistics --name "step05" --input "${run}.vcf.gz" --output "${run}.statistics"

    # Calculate snp call rate and sample call rate (per run)
    java -jar /opt/genomic-utils.jar vcf-quality-control \
        ${run}.vcf.gz \
        --minSnpCallRate ${params.minSnpCallRate}  \
        --minSampleCallRate ${params.minSampleCallRate}  \
        --chunkSize ${params.chunkSize} \
        --output ${run}.qc
    """

}