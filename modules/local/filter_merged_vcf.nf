process FILTER_MERGED_VCF() {

    input:
    path(vcf_file)
    path(vcf_file_index)

    output:
    path("${params.project}.vcf.gz"), emit: vcf_file
    path("${params.project}.vcf.gz.tbi"), emit: vcf_file_index
    path("${params.project}.qc.*"), emit: vcf_file_statistics
    path("${params.project}.statistics"), emit: filter_statistics

    script:
    """

    # Filter by snp call rate and by sample call rate
    vcftools --gzvcf ${vcf_file}  \
        --maf ${params.maf} \
        --hwe ${params.hwe} \
        --recode --stdout | bgzip -c > ${params.project}.filtered.vcf.gz

    java -jar /opt/genomic-utils.jar vcf-statistics --name "maf-hwe"  --input ${params.project}.filtered.vcf.gz --output ${params.project}.statistics

    # Calculate snp call rate and sample call rate
    java -jar /opt/genomic-utils.jar vcf-quality-control \
        ${params.project}.filtered.vcf.gz \
        --minSnpCallRate ${params.minSnpCallRate}  \
        --minSampleCallRate ${params.minSampleCallRate}  \
        --chunkSize ${params.chunkSize} \
        --output ${params.project}.qc

    # Filter by snp call rate and by sample call rate
    vcftools --gzvcf ${params.project}.filtered.vcf.gz  \
        --exclude ${params.project}.qc.snps.excluded  \
        --remove ${params.project}.qc.samples.excluded  \
        --recode --stdout | bgzip -c > ${params.project}.vcf.gz

    tabix ${params.project}.vcf.gz

    java -jar /opt/genomic-utils.jar vcf-statistics --name "final" --input ${params.project}.vcf.gz --output ${params.project}.statistics

    """

}