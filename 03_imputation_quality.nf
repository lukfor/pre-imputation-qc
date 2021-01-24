params.project = "test-gwas"
params.output = "tests/output"
params.chip = "GSAMD-24v3-0-EA_20034606_A1.b37"

params.stepInput = "${params.output}/02_imputation"
params.stepOutput = "${params.output}/03_imputation_quality"

imputation_quality_report = file("$baseDir/reports/03_imputation_quality.Rmd")

InfoFileStatistics = "$baseDir/bin/InfoFileStatistics.java"

// load all laser files
info_viles_ch = Channel.fromPath("${params.stepInput}/*.info.gz")

process calcRsqMean {

  publishDir "$params.stepOutput", mode: 'copy'

  input:
    file info_files from info_viles_ch.collect()

  output:
    file "${params.project}.qualities.txt" into quality_file_ch

  """
  # Calculate mean_Rsq based on MAF
  jbang ${InfoFileStatistics} \
    ${info_files} \
    --output ${params.project}.qualities.txt
  """

}

process createReport {

  publishDir "$params.output", mode: 'copy'

  input:
    file quality_file from quality_file_ch.collect()

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require( 'rmarkdown' ); render('${imputation_quality_report}',
   params = list(
     project = '${params.project}',
     qualities = '${quality_file}'
   ),
   knit_root_dir='\$PWD', output_file='\$PWD/03_imputation_quality.html')"
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
