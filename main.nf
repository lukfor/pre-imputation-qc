#!/usr/bin/env nextflow

include { CLEAN_SAMPLE_IDS } from "./modules/local/clean_sample_ids"
include { EXCLUDE_SAMPLES } from "./modules/local/exclude_samples"
include { FILTER_AND_FIX_STRAND_FLIPS } from "./modules/local/filter_and_fix_strand_flips"
include { MERGE_VCF_FILES } from "./modules/local/merge_vcf_files"
include { FILTER_MERGED_VCF } from "./modules/local/filter_merged_vcf"
include { CREATE_FINAL_PLINK } from "./modules/local/create_final_plink"
include { SPLIT_INTO_CHROMOSOMES } from "./modules/local/split_into_chromosomes"
include { CREATE_REPORT } from "./modules/local/create_report"


requiredParams = [
    'project', 'input',
    'output', 'chip',
    'build', 'strand_file',
    'refalt_file'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}


// load all plink files from input folder
plink_files = Channel.fromFilePairs(
    params.input,
    size: 2,
    flat: true,
    checkIfExists: true
)
chromosomes = Channel.of(1..22)


workflow {

    if (params.cleanSampleIds) {
        CLEAN_SAMPLE_IDS (
            plink_files
        ) 
        plink_files = CLEAN_SAMPLE_IDS.out.vcf_file
    }

    if (params.excludeSamples != null) {
        EXCLUDE_SAMPLES (
            plink_files,
            file(params.excludeSamples, checkIfExists: true)
        )
        plink_files = EXCLUDE_SAMPLES.out.vcf_file
    }

    FILTER_AND_FIX_STRAND_FLIPS (
        plink_files,
        file(params.strand_file, checkIfExists: true),
        file(params.refalt_file, checkIfExists: true)
    )

    MERGE_VCF_FILES (
        FILTER_AND_FIX_STRAND_FLIPS.out.vcf_files.collect(),
        FILTER_AND_FIX_STRAND_FLIPS.out.vcf_files_index.collect()
    )

    FILTER_MERGED_VCF (
        MERGE_VCF_FILES.out.vcf_file.collect(),
        MERGE_VCF_FILES.out.vcf_file_index.collect()
    )

    CREATE_FINAL_PLINK (
        FILTER_MERGED_VCF.out.vcf_file
    )

    SPLIT_INTO_CHROMOSOMES (
        chromosomes,
        FILTER_MERGED_VCF.out.vcf_file.collect(),
        FILTER_MERGED_VCF.out.vcf_file_index.collect()
    )

    CREATE_REPORT(
        MERGE_VCF_FILES.out.vcf_file_statistics.collect(),
        FILTER_MERGED_VCF.out.vcf_file_statistics.collect(),
        FILTER_AND_FIX_STRAND_FLIPS.out.samples_runs.collect(),
        FILTER_AND_FIX_STRAND_FLIPS.out.snps_runs.collect(),
        FILTER_AND_FIX_STRAND_FLIPS.out.filter_statistics.collect(),
        FILTER_MERGED_VCF.out.filter_statistics.collect()
    )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
