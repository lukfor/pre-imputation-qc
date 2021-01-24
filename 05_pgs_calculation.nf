params.project = "test-gwas"
params.output = "tests/output"
params.chip = "GSAMD-24v3-0-EA_20034606_A1.b37"
params.pgs_scores = "PGS000115,PGS000013,PGS000340,PGS000015,PGS000004"

params.stepInput = "${params.output}/02_imputation"
params.stepOutput = "${params.output}/05_pgs_calculation"


imputed_vcf_files_ch = Channel.fromPath("${params.stepInput}/*.dose.vcf.gz")

process calcScores {

  publishDir "$params.stepOutput", mode: 'copy'

  input:
    file imputed_vcf_files from imputed_vcf_files_ch.collect()

  output:
    file "${params.project}.scores.txt" into results_ch
    file "${params.project}.scores.html" into report_ch

  """
  # Download latest meta data
  wget https://www.pgscatalog.org/rest/score/all -O pgs-catalog.json

  # TODO: use samples file from outliers

  pgs-calc ${imputed_vcf_files} \
    --ref ${params.pgs_scores} \
    --out ${params.project}.scores.txt \
    --report-html ${params.project}.scores.html \
    --meta pgs-catalog.json \
    --no-ansi

  """

}

process createReport {

  publishDir "$params.output", mode: 'copy'

  input:
    file report from report_ch.collect()

  output:
    file "*.html" into final_report_ch

  """
  cp ${report} 05_pgs_calculation.html
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
