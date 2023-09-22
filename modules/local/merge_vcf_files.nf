process MERGE_VCF_FILES {

    input:
    path(vcf_files)
    path(vcf_files_index)

    output:
    path("${params.project}.merged.vcf.gz"), emit: vcf_file
    path("${params.project}.merged.vcf.gz.tbi"), emit: vcf_file_index
    path("${params.project}.merged.statistics"), emit: vcf_file_statistics

    script:
    """

    # if contains a spaceh --> multiple files --> merge needed
    if [[ "${vcf_files}" = *" "* ]]; then
        # TODO: -m id still needed with extract? use -O z to avoid bgzip.
        bcftools merge -m id ${vcf_files} -O v > ${params.project}.merged.vcf
        bgzip ${params.project}.merged.vcf
    else
        cp ${vcf_files} ${params.project}.merged.vcf.gz
    fi

    tabix ${params.project}.merged.vcf.gz

    java -jar /opt/genomic-utils.jar vcf-statistics --name "merged" --input ${params.project}.merged.vcf.gz --output ${params.project}.merged.statistics

    """

}