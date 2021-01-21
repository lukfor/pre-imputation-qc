params.project = "test-gwas"
params.output = "tests/output"
params.chip = "GSAMD-24v3-0-EA_20034606_A1.b37"
params.pca_reference_samples = "$baseDir/data/hgdp.txt"
params.pca_reference_population = "Europe"
params.pca_max_pc = 3
params.pca_max_sd = 6

params.stepInput = "${params.output}/02_imputation"
params.stepOutput = "${params.output}/03_post_imputation"

post_imputation_report = file("$baseDir/reports/03_post_imputation.Rmd")

FilterPCA = "$baseDir/bin/FilterPCA.java"

// load all laser files
laser_files_ch = Channel.fromPath("${params.stepInput}/*.txt")
laser_files2_ch = Channel.fromPath("${params.stepInput}/*.txt")
pca_reference_samples_ch = file(params.pca_reference_samples)

process detectOutliers {

  publishDir "$params.stepOutput", mode: 'copy'

  input:
    file laser_files from laser_files_ch.collect()
    file pca_reference_samples from pca_reference_samples_ch

  output:
    file "${params.project}*.txt" into outlier_files_ch

  """
  # Calculate pca clusters
  jbang ${FilterPCA} \
    --reference-pc reference_pc.txt \
    --reference-samples ${pca_reference_samples} \
    --population ${params.pca_reference_population} \
    --max-pc ${params.pca_max_pc} \
    --max-sd ${params.pca_max_sd} \
    --study-pc study_pc.txt \
    --output ${params.project}
  """

}

process createReport {

  publishDir "$params.output", mode: 'copy'

  input:
    file laser_files from laser_files2_ch.collect()
    file outlier_files from outlier_files_ch.collect()
    file pca_reference_samples from pca_reference_samples_ch

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require( 'rmarkdown' ); render('${post_imputation_report}',
   params = list(
     project = '${params.project}',
     reference_pc = 'reference_pc.txt',
     reference_pc_var = 'reference_pc_var.txt',
     reference_samples = '${pca_reference_samples}',
     study_pc = 'study_pc.txt',
     pca_reference_population = '${params.pca_reference_population}',
     pca_max_pc = '${params.pca_max_pc}',
     pca_max_sd = '${params.pca_max_sd}',
     study_outliers = '${params.project}.outliers.txt',
     reference_cluster = '${params.project}.clusters.txt'
   ),
   knit_root_dir='\$PWD', output_file='\$PWD/03_post_imputation.html')"
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
