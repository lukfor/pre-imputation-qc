params.stepInput = "${params.output}/typed/vcf/*chr*.vcf.gz"
params.stepOutput = "${params.output}/imputed/${params.refPanel}"

imputation_quality_report = file("$baseDir/reports/03_imputation_quality.Rmd")

InfoFileStatistics = "$baseDir/bin/InfoFileStatistics.java"


// load all vcf files from input folder
vcf_files_ch = Channel.fromPath("${params.stepInput}")

if (params.token == "") {
   exit 1, "Parameter 'token' is required"
}

process imputeGenotypes {

  publishDir "${params.stepOutput}/vcf", mode: 'copy',
    saveAs: { "${file(it).getName()}" }

  input:
    file vcf_files from vcf_files_ch.collect()

  output:
    file "**/local/*.dose.vcf.gz" into imputed_vcf_files_ch
    file "**/local/*.info.gz" into info_files_ch
    file "**/qcreport/*.html"
    file "**/statisticDir/*.txt"

  """
  # configure imputationbot
  echo -e "-  hostname: ${params.imputation_server}\n   token: ${params.imputation_token}\n" > ~/.imputationbot/imputationbot.instances

  imputationbot impute \
    --files ${vcf_files} \
    --refpanel ${params.imputation_reference_panel} \
    --population ${params.imputation_population} \
    --autoDownload \
    --password ${params.imputation_password} \
    --build ${params.build}

  #TODO: delete zip files

  """

}

process calcRsqMean {

  input:
    file info_files from info_files_ch.collect()

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

  publishDir "${params.stepOutput}", mode: 'copy'

  input:
    file quality_file from quality_file_ch.collect()
    file imputation_quality_report

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require( 'rmarkdown' ); render('${imputation_quality_report}',
   params = list(
     project = '${params.project}',
     chip = '${params.chip}',
     reference_panel = '${params.imputation_reference_panel}'
     qualities = '${quality_file}'
   ),
   knit_root_dir='\$PWD', output_file='\$PWD/imputation_quality.html')"
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
