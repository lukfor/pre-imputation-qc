params.stepInput = "${params.input}"
params.stepOutput = "${params.output}/genotyped"

VcfQualityControl = "$baseDir/bin/VcfQualityControl.java"
VcfStatistics = "$baseDir/bin/VcfStatistics.java"
pre_imputation_report = file("$baseDir/reports/01_pre_imputation.Rmd")

pca_report = file("$baseDir/reports/04_pca_smartpca.Rmd")
ibd_report = file("$baseDir/reports/05_ibd.Rmd")
high_ld_file = file("$baseDir/data/high-ld.txt")


requiredParams = [
    'project', 'input',
    'output', 'chip',
    'build'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}


// TODO: check strand/refalt file

// load all plink files from input folder
Channel.fromFilePairs("${params.stepInput}").set {plink_files_ch}
strand_file_ch = file(params.strand_file)
refalt_file_ch = file(params.refalt_file)
chromosomes_ch = Channel.of(1..22)


if (params.cleanSampleIds) {

  process cleanSampleIds {

    input:
      set filename, file(map_file) from plink_files_ch

    output:
      tuple val("${filename}.cleaned"), file("${filename}.cleaned.*") into plink_files_cleaned_ch

    """
    # Step 0: replace all spaces with underscore (e.g. spaces in Sample IDs)
    sed -e 's/ /_/g' ${filename}.ped > ${filename}.cleaned.ped
    cp ${filename}.map ${filename}.cleaned.map
    """

  }

} else {

  plink_files_cleaned_ch = plink_files_ch

}

if (params.excludeSamples != null) {

  exclude_samples_file = file(params.excludeSamples)

  process excludeSamples {

    input:
      set filename, file(map_file) from plink_files_cleaned_ch
      file exclude_samples_file

    output:
      tuple val("${filename}.excluded"), file("${filename}.excluded.*") into plink_files_excluded_ch

    """
    plink --file ${filename} \
      --remove ${exclude_samples_file} \
      --recode \
      --tab \
      --out ${filename}.excluded
    """

  }

} else {

  plink_files_excluded_ch = plink_files_cleaned_ch

}


process filterAndFixStrandFlips {

  input:
    set filename, file(map_file) from plink_files_excluded_ch
    file strand_file from strand_file_ch
    file refalt_file from refalt_file_ch

  output:
    file "${filename}.vcf.gz" into vcf_files_ch
    file "${filename}.vcf.gz.tbi" into vcf_files_index_ch
    file "${filename}.qc.samples" into samples_runs_ch
    file "${filename}.qc.snps" into snps_runs_ch
    file "*.statistics" into filter_statistics_ch

  """

  # count all lines in map file and ped file as write as step0 to .statistics
  total_samples=\$(cat ${filename}.ped | wc -l)
  total_snps=\$(cat ${filename}.map | wc -l)
  printf "name samples snps\n" > ${filename}.statistics
  printf "step00 \${total_samples} \${total_snps}" >> ${filename}.statistics

  # Step 1: Remove all indels, "I" and "D"
  plink --file ${filename} \
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
     --a2-allele ${refalt_file} \
     --recode ${params.useDoubleId ? 'vcf-iid' : 'vcf'} \
     --out ${filename}.harmonized

  # remove all variants that have no ref allele inside

  grep  "Warning: Impossible A2 allele assignment for variant *" ${filename}.harmonized.log | awk '{print substr(\$NF,1,length(\$NF)-1)}' > ${filename}.harmonized.snps

  vcftools --vcf ${filename}.harmonized.vcf \
    --exclude ${filename}.harmonized.snps \
    --recode --stdout | bgzip -c > ${filename}.vcf.gz

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

  input:
    file vcf_files from vcf_files_ch.collect()
    file vcf_files_index from vcf_files_index_ch.collect()

  output:
    file "${params.project}.merged.vcf.gz" into merged_vcf_file_ch
    file "${params.project}.merged.vcf.gz.tbi" into merged_vcf_file_index_ch
    file "${params.project}.merged.statistics" into merged_vcf_file_statistics

  """

  # if contains a spaceh --> multiple files --> merge needed
  if [[ "${vcf_files}" = *" "* ]]; then
    # TODO: -m id still needed with extract? use -O z to avoid bgzip.
    bcftools merge -m id ${vcf_files} -O v > ${params.project}.merged.vcf
    bgzip ${params.project}.merged.vcf
  else
    cp ${vcf_files} ${params.project}.merged.vcf.gz
  fi

  tabix ${params.project}.merged.vcf.gz

  vcf-statistics "merged" ${params.project}.merged.vcf.gz ${params.project}.merged.statistics

  """

}

process filterMergedVcf() {

  input:
    file merged_vcf_file from merged_vcf_file_ch.collect()
    file merged_vcf_file_index from merged_vcf_file_index_ch.collect()

  output:
    file "${params.project}.vcf.gz" into final_vcf_file_ch
    file "${params.project}.vcf.gz.tbi" into final_vcf_file_index_ch
    file "${params.project}.qc.*" into final_vcf_file_statistics
    file "${params.project}.statistics" into merged_filter_statistics_ch
  """

  # Filter by snp call rate and by sample call rate
  vcftools --gzvcf ${merged_vcf_file}  \
    --maf ${params.maf} \
    --hwe ${params.hwe} \
    --recode --stdout | bgzip -c > ${params.project}.filtered.vcf.gz

  vcf-statistics "maf-hwe" ${params.project}.filtered.vcf.gz ${params.project}.statistics

  # Calculate snp call rate and sample call rate
  jbang ${VcfQualityControl} ${params.project}.filtered.vcf.gz \
    --minSnpCallRate ${params.minSnpCallRate}  \
    --minSampleCallRate ${params.minSampleCallRate}  \
    --chunkSize ${params.chunkSize} \
    --output ${params.project}.qc

  # Filter by snp call rate and by sample call rate
  vcftools --gzvcf ${params.project}.filtered.vcf.gz  \
    --exclude ${params.project}.qc.snps.excluded  \
    --remove ${params.project}.qc.samples.excluded  \
    --recode --stdout | bgzip -c > ${params.project}.vcf.gz

  tabix ${params.project}.vcf.gz

  vcf-statistics "final" ${params.project}.vcf.gz ${params.project}.statistics

  """

}

process createFinalPlink() {

  publishDir "$params.stepOutput/plink", mode: 'copy'

  input:
    file merged_vcf_file from final_vcf_file_ch

  output:
    file "${params.project}.{bim,bed,fam}" into final_plink_file_ch
  """

  plink --vcf ${merged_vcf_file} --double-id --out ${params.project}

  """

}

process splitIntoChromosomes {

  publishDir "$params.stepOutput/vcf", mode: 'copy'

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

  publishDir "$params.stepOutput", mode: 'copy'

  input:
    file stats from merged_vcf_file_statistics.collect()
    file stats2 from final_vcf_file_statistics.collect()
    file samples_runs from samples_runs_ch.collect()
    file snps_runs from snps_runs_ch.collect()
    file filter_statistics from filter_statistics_ch.collect()
    file merged_filter_statistics from merged_filter_statistics_ch.collect()
    file pre_imputation_report

  output:
    file "*.html" into report_ch

  """
  Rscript -e "require( 'rmarkdown' ); render('${pre_imputation_report}',
    params = list(
      project = '${params.project}',
      chip = '${params.chip}',
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
      samples_merged = '${params.project}.merged.statistics'
    ), knit_root_dir='\$PWD', output_file='\$PWD/pre_imputation.html')"
  """

}

if (params.pca_enabled){

  process smartpca {

    input:
      file plink_file from final_plink_file_ch.collect()
      file high_ld_file

    output:
      file "*.{evec,par,out,snps.weights}" into smartpca_files_ch
      file '*.prune.in'
      file '*.set'
      file 'plink.genome' into ibd_estimation_ch

    """
    # Prune, filter vcf and convert to plink (TODO: Check prune parameters and move into own process)
    # We will conduct principle component analysis on genetic variants that are pruned for variants in linkage
    # disequilibrium (LD) with an r2 > 0.2 in a 50kb window

    plink --bfile ${params.project} --double-id --maf 0.01 --hwe 1E-6 --indep-pairwise 50 5 0.2 --out ${params.project}
    plink --bfile ${params.project} --extract ${params.project}.prune.in --double-id --make-bed --out ${params.project}.pruned

    # There are regions of long-range, high linkage diequilibrium in the human genome. These regions should be excluded when performing certain analyses such as principal component analysis on genotype data.

    plink --bfile ${params.project}.pruned --make-set ${high_ld_file} --write-set --out hild
    plink --bfile ${params.project}.pruned --exclude hild.set --make-bed --out ${params.project}.pruned2

    # problem with long id names in bim file --> recreate with new id?
    awk '{print \$1,\$1"_"\$3,\$3,\$4,\$5,\$6}' ${params.project}.bim > ${params.project}.pruned2.updated.bim

    # todo: filter relateness pi_hat > 0.1875 (see wuttke et.al )
    plink --bfile ${params.project}.pruned2 --genome


    # TODO: create pedind file from fam file?

    echo "genotypename: ${params.project}.pruned2.bed" > ${params.project}.par
    echo "snpname: ${params.project}.pruned2.updated.bim" >> ${params.project}.par
    echo "indivname: ${params.project}.pruned2.fam" >> ${params.project}.par
    echo "evecoutname: ${params.project}.evec" >> ${params.project}.par
    echo "evaloutname: ${params.project}.eval" >> ${params.project}.par
    echo "snpweightoutname: ${params.project}.snps.weights" >> ${params.project}.par
    echo "altnormstyle: NO" >> ${params.project}.par
    echo "numoutlieriter: 0" >> ${params.project}.par
    echo "familynames: NO" >> ${params.project}.par
    echo "numoutevec: ${params.pca_max_pc }" >> ${params.project}.par

    smartpca -p ${params.project}.par >  ${params.project}.out
    """

  }

  process createRelatenessReport {

    publishDir "${params.stepOutput}/ibd", mode: 'copy'

    input:
      file ibd_estimation_file from ibd_estimation_ch
      file ibd_report

    output:
      file "plink.genome"
      file "*.html"

    """
    Rscript -e "require( 'rmarkdown' ); render('${ibd_report}',
      params = list(
        genome_filename = '${ibd_estimation_file}'
      ), knit_root_dir='\$PWD', output_file='\$PWD/ibd_smartpca.html')"
    """

  }

  process createPcaReport {

    publishDir "${params.stepOutput}/pca", mode: 'copy'

    input:
      file stats from smartpca_files_ch.collect()
      file pca_report

    output:
      file "${params.project}.pca.txt"
      file "*.html"

    """
    Rscript -e "require( 'rmarkdown' ); render('${pca_report}',
      params = list(
        evec_filename = '${params.project}.evec',
        output_filename = '${params.project}.pca.txt',
        snps_weights_filename = '${params.project}.snps.weights',
        max_pc = '${params.pca_max_pc}'
      ), knit_root_dir='\$PWD', output_file='\$PWD/pca_smartpca.html')"
    """

  }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
