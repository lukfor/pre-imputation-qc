process CREATE_FINAL_PLINK() {

    publishDir "$params.output/plink", mode: 'copy'

    input:
    path(vcf_file)

    output:
    file "${params.project}.{bim,bed,fam}"

    script:
    """
    plink --vcf ${vcf_file} --double-id --out ${params.project}
    """

}