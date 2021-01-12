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
  # Step 1: replace all spaces with underscore (e.g. spaces in Sample IDs)
  sed -e 's/ /_/g' ${filename}.ped > ${filename}.step01.ped
  cp ${filename}.map ${filename}.step01.map

  # TODO: count all lines in map file


  # Step 2: Remove all indels, "I" and "D"

  # Write list of all indel calls (D/I allele codes)
  plink --file ${filename}.step01 --list-23-indels --out ${filename}.step01
  plink --file ${filename}.step01 \
    --exclude ${filename}.step01.indel \
    --make-bed \
    --out ${filename}.step02


  #plink --file ${filename}.step01 \
#    --snps-only just-acgt \
    #--make-bed \
    #--out ${filename}.step02

  plink-statistics "step02_remove_all_indels" ${filename}.step02 ${filename}.statistics


  # Step 3: Remove all non-autosomale SNPs
  plink --bfile ${filename}.step02 \
    --chr 1-22 \
    --make-bed \
    --out ${filename}.step03

  plink-statistics "step03_remove_non_autosomale" ${filename}.step03 ${filename}.statistics


  # Step 4: Update strand flips (https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh)
  update_build.sh \
    ${filename}.step03 \
    ${strand_file} \
    ${filename}.step04

  plink-statistics "step04_update_strand_flips" ${filename}.step04 ${filename}.statistics

  # Step 5: Remvoe all autosomale snps after update strand flips
  plink --bfile ${filename}.step04 \
    --chr 1-22 \
    --make-bed \
    --out ${filename}.step05

    plink-statistics "step05_remove_non_autosomale" ${filename}.step05 ${filename}.statistics


  # Step 6: Harmonize ref/alt alleles and retain only SNPs in the refalt file.
  plink --bfile ${filename}.step05 \
    --extract ${refalt_file} \
    --reference-allele ${refalt_file} \
    --recode vcf \
    --out ${filename}

  bgzip ${filename}.vcf
  tabix ${filename}.vcf.gz

  vcf-statistics "step06_harmonize_ref_alt_alleles" ${filename}.vcf.gz ${filename}.statistics


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
    file "${params.project}.vcf.gz" into merged_vcf_file_ch
    file "${params.project}.vcf.gz.tbi" into merged_vcf_file_index_ch
    file "${params.project}.qc.*" into merged_vcf_statistics
    file "${params.project}.statistics" into merged_vcf_statistics2

  """
  # TODO: check length of vc_files. if only one element, no merge needed, copy only.

  # TODO: -m id still needed with extract? use -O z to avoid bgzip.
  bcftools merge -m id ${vcf_files} -O v > ${params.project}.unfiltered.vcf
  bgzip ${params.project}.unfiltered.vcf
  tabix ${params.project}.unfiltered.vcf.gz

  vcf-statistics "unfiltered" ${params.project}.unfiltered.vcf.gz ${params.project}.statistics

  # Calculate snp call rate and sample call rate
  jbang ${VcfQualityControl} ${params.project}.unfiltered.vcf.gz \
    --minSnpCallRate ${params.minSnpCallRate}  \
    --minSampleCallRate ${params.minSampleCallRate}  \
    --chunkSize ${params.chunkSize} \
    --output ${params.project}.qc

  # Filter by snp call rate and by sample call rate
  vcftools --gzvcf ${params.project}.unfiltered.vcf.gz  \
    --exclude ${params.project}.qc.snps.excluded  \
    --remove ${params.project}.qc.samples.excluded  \
    --recode --stdout | bgzip -c > ${params.project}.vcf.gz

  tabix ${params.project}.vcf.gz

  vcf-statistics "filtered" ${params.project}.vcf.gz ${params.project}.statistics

  """

}


process splitIntoChromosomes {

  publishDir "$params.output", mode: 'copy'

  input:
    val chromosome from chromosomes_ch
    file vcf_file from merged_vcf_file_ch.collect()
    file vcf_file_index_index from merged_vcf_file_index_ch.collect()

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
    file stats from merged_vcf_statistics.collect()
    file stats2 from merged_vcf_statistics2.collect()
    file samples_runs from samples_runs_ch.collect()
    file snps_runs from snps_runs_ch.collect()
    file filter_statistics from filter_statistics_ch.collect()
    file pre_imputation_report

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require( 'rmarkdown' ); render('${pre_imputation_report}', params = list(project = '${params.project}', chip = '${params.chip}', samples = '${samples_runs}', snps = '${snps_runs}', filter_statistics = '${filter_statistics}', samples_excluded = '${params.project}.qc.samples.excluded', samples_merged = '${params.project}.qc.samples'), knit_root_dir='\$PWD', output_file='\$PWD/pre-imputation-report.html')"
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
