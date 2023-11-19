process CREATE_REPORT {

  publishDir "$params.output", mode: 'copy'

  input:
    path(stats)
    path(stats2)
    path(samples_runs)
    path(snps_runs)
    path(filter_statistics)
    path(merged_filter_statistics)

  output:
    file "*.html"

  """
  Rscript - <<EOF
    require( 'rmarkdown' ); render('`which report.Rmd`',   
      params = list(
        project = '${params.project}',
        date = '${params.project_date}',
        chip = '${params.chip}',
        build = '${params.build}',
        maf = '${params.maf}',
        hwe = '${params.hwe}',
        snp_call_rate = '${params.minSnpCallRate}',
        sample_call_rate = '${params.minSampleCallRate}',
        samples = '${samples_runs}',
        snps = '${snps_runs}',
        filter_statistics = '${filter_statistics}',
        merged_filter_statistics = '${merged_filter_statistics}',
        samples_excluded = '${params.project}.qc.samples.excluded',
        snps_excluded = '${params.project}.qc.snps.excluded',
        samples_final = '${params.project}.qc.samples',
        snps_final = '${params.project}.qc.snps',
        samples_merged = '${params.project}.merged.statistics',
        version  ='${workflow.manifest.version}'
      ), knit_root_dir='\$PWD', output_file='\$PWD/report.html')
  EOF
  """

}