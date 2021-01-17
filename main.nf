params.project = "test-gwas"
params.input = "tests/input/*/*.{map,ped}"
params.output = "tests/output"
params.chip = "GSAMD-24v3-0-EA_20034606_A1.b37"

params.chunkSize= 20000000
params.minSampleCallRate = 0.5
params.minSnpCallRate = 0.9

params.strand_file = "$baseDir/data/${params.chip}.strand"
params.refalt_file = "$baseDir/data/${params.chip}.RefAlt"


VcfQualityControl = "$baseDir/bin/VcfQualityControl.java"
VcfStatistics = "$baseDir/bin/VcfStatistics.java"
pre_imputation_report = file("$baseDir/reports/pre-imputation.Rmd")

if (!params.input) {
    exit 1, "Plink files not specified"
}

// TODO: check strand/refalt file

// load all plink files from input folder
Channel.fromFilePairs("${params.input}").set {plink_files_ch}
strand_file_ch = file(params.strand_file)
refalt_file_ch = file(params.refalt_file)
chromosomes_ch = Channel.of(1..22)


process filterAndFixStrandFlips {

  publishDir "$params.output/single", mode: 'copy'

  input:
    set filename, file(map_file) from plink_files_ch
    file strand_file from strand_file_ch
    file refalt_file from refalt_file_ch

  output:
    file "${filename}.vcf.gz" into vcf_files_ch
    file "${filename}.vcf.gz.tbi" into vcf_files_index_ch
    file "${filename}.qc.samples" into samples_runs_ch
    file "${filename}.qc.snps" into snps_runs_ch
    file "*.statistics" into filter_statistics_ch

  """
  # Step 0: replace all spaces with underscore (e.g. spaces in Sample IDs)
  sed -e 's/ /_/g' ${filename}.ped > ${filename}.step00.ped
  cp ${filename}.map ${filename}.step00.map

  # count all lines in map file and ped file as write as step0 to .statistics
  total_samples=\$(cat ${filename}.step00.ped | wc -l)
  total_snps=\$(cat ${filename}.step00.map | wc -l)
  printf "name samples snps\n" > ${filename}.statistics
  printf "step00 \${total_samples} \${total_snps}" >> ${filename}.statistics

  # Step 1: Remove all indels, "I" and "D"
  plink --file ${filename}.step00 \
    --snps-only just-acgt \
    --make-bed \
    --out ${filename}.step01

  plink-statistics "step01" ${filename}.step01 ${filename}.statistics


  # Step 2: Remove all non-autosomale SNPs
  plink --bfile ${filename}.step01 \
    --chr 1-22 \
    --make-bed \
    --out ${filename}.step02

  plink-statistics "step02" ${filename}.step02 ${filename}.statistics


  # Step 3: Update strand flips (https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh)
  update_build.sh \
    ${filename}.step02 \
    ${strand_file} \
    ${filename}.step03

  plink-statistics "step03" ${filename}.step03 ${filename}.statistics

  # Step 4: Remove all autosomale snps after update strand flips
  plink --bfile ${filename}.step03 \
    --chr 1-22 \
    --make-bed \
    --out ${filename}.step04

    plink-statistics "step04" ${filename}.step04 ${filename}.statistics


  # Step 5: Harmonize ref/alt alleles and retain only SNPs in the refalt file.
  plink --bfile ${filename}.step04 \
    --extract ${refalt_file} \
    --reference-allele ${refalt_file} \
    --recode vcf \
    --out ${filename}

  bgzip ${filename}.vcf
  tabix ${filename}.vcf.gz

  vcf-statistics "step05" ${filename}.vcf.gz ${filename}.statistics


  # Calculate snp call rate and sample call rate (per run)
  jbang ${VcfQualityControl} ${filename}.vcf.gz \
    --minSnpCallRate ${params.minSnpCallRate}  \
    --minSampleCallRate ${params.minSampleCallRate}  \
    --chunkSize ${params.chunkSize} \
    --output ${filename}.qc
  """

}


process mergeVcfFiles() {

  publishDir "$params.output", mode: 'copy'

  input:
    file vcf_files from vcf_files_ch.collect()
    file vcf_files_index from vcf_files_index_ch.collect()

  output:
    file "${params.project}.merged.vcf.gz" into merged_vcf_file_ch
    file "${params.project}.merged.vcf.gz.tbi" into merged_vcf_file_index_ch
    file "${params.project}.merged.statistics" into merged_vcf_file_statistics

  """

  # TODO: check length of vc_files. if only one element, no merge needed, copy only.

  # TODO: -m id still needed with extract? use -O z to avoid bgzip.
  bcftools merge -m id ${vcf_files} -O v > ${params.project}.merged.vcf
  bgzip ${params.project}.merged.vcf
  tabix ${params.project}.merged.vcf.gz

  vcf-statistics "merged" ${params.project}.merged.vcf.gz ${params.project}.merged.statistics

  """

}

process filterMergedVcf() {

  publishDir "$params.output", mode: 'copy'

  input:
    file merged_vcf_file from merged_vcf_file_ch.collect()
    file merged_vcf_file_index from merged_vcf_file_index_ch.collect()

  output:
    file "${params.project}.vcf.gz" into final_vcf_file_ch
    file "${params.project}.vcf.gz.tbi" into final_vcf_file_index_ch
    file "${params.project}.qc.*" into final_vcf_file_statistics

  """

  # Calculate snp call rate and sample call rate
  jbang ${VcfQualityControl} ${merged_vcf_file} \
    --minSnpCallRate ${params.minSnpCallRate}  \
    --minSampleCallRate ${params.minSampleCallRate}  \
    --chunkSize ${params.chunkSize} \
    --output ${params.project}.qc

  # Filter by snp call rate and by sample call rate
  vcftools --gzvcf ${merged_vcf_file}  \
    --exclude ${params.project}.qc.snps.excluded  \
    --remove ${params.project}.qc.samples.excluded  \
    --recode --stdout | bgzip -c > ${params.project}.vcf.gz

  tabix ${params.project}.vcf.gz

  vcf-statistics "final" ${params.project}.vcf.gz ${params.project}.statistics

  """

}

process splitIntoChromosomes {

  publishDir "$params.output", mode: 'copy'

  input:
    val chromosome from chromosomes_ch
    file vcf_file from final_vcf_file_ch.collect()
    file vcf_file_index_index from final_vcf_file_index_ch.collect()

  output:
    file "${params.project}.chr*.vcf.gz" into chr_vcf_file_ch

  """
  # TODO: use -O z to avoid bgzip?
  bcftools view -r ${chromosome} ${vcf_file} > ${params.project}.chr${chromosome}.vcf
  bgzip ${params.project}.chr${chromosome}.vcf
  tabix ${params.project}.chr${chromosome}.vcf.gz
  """

}


process createReport {

  publishDir "$params.output", mode: 'copy'

  input:
    file stats from merged_vcf_file_statistics.collect()
    file stats2 from final_vcf_file_statistics.collect()
    file samples_runs from samples_runs_ch.collect()
    file snps_runs from snps_runs_ch.collect()
    file filter_statistics from filter_statistics_ch.collect()
    file pre_imputation_report

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require( 'rmarkdown' ); render('${pre_imputation_report}', params = list(project = '${params.project}', chip = '${params.chip}', samples = '${samples_runs}', snps = '${snps_runs}', filter_statistics = '${filter_statistics}', samples_excluded = '${params.project}.qc.samples.excluded', snps_excluded = '${params.project}.qc.snps.excluded', samples_final = '${params.project}.qc.samples', snps_final = '${params.project}.qc.snps', samples_merged = '${params.project}.merged.statistics'), knit_root_dir='\$PWD', output_file='\$PWD/pre-imputation-report.html')"
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
