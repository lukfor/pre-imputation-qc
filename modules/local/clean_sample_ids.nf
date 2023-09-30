process CLEAN_SAMPLE_IDS {

    input:
    tuple val(run), val(filename), path(map_file), path(ped_file)

    output:
    tuple val(run), val("${run}.cleaned"), file("${run}.cleaned.map"), file("${run}.cleaned.ped"), emit: vcf_file

    script:
    """
    # Step 0: replace all spaces with underscore (e.g. spaces in Sample IDs)
    sed -e 's/ /_/g' "${filename}.ped" > "${run}.cleaned.ped"
    cp "${filename}.map" "${run}.cleaned.map"
    """

}