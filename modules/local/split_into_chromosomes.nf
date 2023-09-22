process SPLIT_INTO_CHROMOSOMES {

    publishDir "$params.output/vcf", mode: 'copy'

    input:
    val(chromosome)
    path(vcf_file)
    path(vcf_file_index)

    output:
    path "${params.project}.chr*.vcf.gz"

    script:
    """
    # TODO: use -O z to avoid bgzip?
    bcftools view -r ${chromosome} ${vcf_file} > ${params.project}.chr${chromosome}.vcf
    bgzip ${params.project}.chr${chromosome}.vcf
    tabix ${params.project}.chr${chromosome}.vcf.gz
    """

}
