params.project = "test-gwas"
params.input = "tests/input/*/*.{map,ped}"
params.output = "tests/output"
params.chip = "GSAMD-24v3-0-EA_20034606_A1.b37"

params.chunkSize= 20000000
params.minSampleCallRate = 0.5
params.minSnpCallRate = 0.9

params.strand_file = "data/${params.chip}.strand"
params.refalt_file = "data/${params.chip}.RefAlt"


VcfQualityControl = "$baseDir/bin/VcfQualityControl.java"

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

  """
  # replace all spaces with underscore (e.g. spaces in Sample IDs)
  sed -e 's/ /_/g' ${filename}.ped > ${filename}.fixed.ped
  cp ${filename}.map ${filename}.fixed.map

  # TODO: write statistics before filtering. convert immideatly to vcf and update "update_build.sh" script?

  # Remove all indels, "I" and "D"
  plink --file ${filename}.fixed \
    --snps-only just-acgt \
    --make-bed \
    --out ${filename}.binary

  # Remove all non-autosomale SNPs
  plink --bfile ${filename}.binary \
    --chr 1-22 \
    --make-bed \
    --out ${filename}.autosomes

  # https://www.well.ox.ac.uk/~wrayner/strand/update_build.sh
  update_build.sh \
    ${filename}.autosomes \
    ${strand_file} \
    ${filename}.autosomes.strand

  # Harmonize ref/alt alleles and retain only SNPs in the refalt file
  plink --bfile ${filename}.autosomes.strand \
    --chr 1-22 \
    --extract ${refalt_file} \
    --reference-allele ${refalt_file} \
    --recode vcf \
    --out ${filename}

  bgzip ${filename}.vcf
  tabix ${filename}.vcf.gz

  # TODO: write statistics after filtering

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

  """
  # TODO: check length of vc_files. if only one element, no merge needed, copy only.

  # TODO: -m id still needed with extract? use -O z to avoid bgzip.
  bcftools merge -m id ${vcf_files} -O v > ${params.project}.unfiltered.vcf
  bgzip ${params.project}.unfiltered.vcf
  tabix ${params.project}.unfiltered.vcf.gz

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
    file samples_runs from samples_runs_ch.collect()
    file snps_runs from snps_runs_ch.collect()

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require( 'rmarkdown' ); render('$baseDir/reports/pre-imputation.Rmd', params = list(project = '${params.project}', chip = '${params.chip}', samples = '${samples_runs}', snps = '${snps_runs}', samples_excluded = '${params.project}.qc.samples.excluded'), knit_root_dir='\$PWD', output_file='\$PWD/pre-imputation-report.html')"
  """

}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
