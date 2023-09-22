process EXCLUDE_SAMPLES {

    input:
    tuple val(filename), path(map_file), path(ped_file)
    path(exclude_samples_file)

    output:
    tuple val("${filename}.excluded"), file("${filename}.excluded.*"), emit: vcf_file

    script:
    """
    plink --file ${filename} \
      --remove ${exclude_samples_file} \
      --recode \
      --tab \
      --out ${filename}.excluded
    """

}