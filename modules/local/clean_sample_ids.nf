process CLEAN_SAMPLE_IDS {

    input:
    tuple val(filename), path(map_file), path(ped_file)

    output:
    tuple val("${filename}.cleaned"), file("${filename}.cleaned.*"), emit: vcf_file

    script:
    """
    # Step 0: replace all spaces with underscore (e.g. spaces in Sample IDs)
    sed -e 's/ /_/g' ${filename}.ped > ${filename}.cleaned.ped
    cp ${filename}.map ${filename}.cleaned.map
    """

}