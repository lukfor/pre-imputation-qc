process EXCLUDE_SAMPLES {

    input:
    tuple val(run), val(filename), path(map_file), path(ped_file)
    path(exclude_samples_file)

    output:
    tuple val(run), val("${run}.excluded"), file("${run}.excluded.map"), file("${run}.excluded.ped"), emit: vcf_file

    script:
    """
    plink --file "${filename}" \
      --remove ${exclude_samples_file} \
      --recode \
      --tab \
      --out ${run}.excluded
    """

}