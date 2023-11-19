process GET_SAMPLES {
    
    input:
    tuple(path(vcf_file), val(source))

    output:
    path("${vcf_file}.samples")

    script:
    """
    bcftools query -l ${vcf_file} | while read line; do echo "\$line;${source}"; done >> ${vcf_file}.samples
    """

}