params.project = "test-gwas"
params.output = "tests/output"
params.refPanel = "1000g-phase-3-v5"
params.population = "eur"
params.password = "lukas_48318786414"

params.stepInput = "${params.output}/01_pre_imputation"
params.stepOutput = "${params.output}/02_imputation"

// load all plink files from input folder
vcf_files_ch = Channel.fromPath("${params.stepInput}/*chr*.vcf.gz")

process imputeGenotypes {

  publishDir "$params.output", mode: 'copy'

  input:
    file vcf_files from vcf_files_ch.collect()

  output:
    file "**/local/*.dose.vcf.gz" into imputed_vcf_files_ch
    file "**/local/*.info.gz" into info_files_ch
  """
  # TODO: export IMPUTATION_BOT token?
  imputationbot impute \
    --files ${vcf_files} \
    --refpanel ${params.refPanel} \
    --population ${params.population} \
    --autoDownload \
    --password ${params.password}
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
