#!/usr/bin/env nextflow

include { CLEAN_SAMPLE_IDS } from "./modules/local/clean_sample_ids"
include { EXCLUDE_SAMPLES } from "./modules/local/exclude_samples"
include { FILTER_AND_FIX_STRAND_FLIPS } from "./modules/local/filter_and_fix_strand_flips"
include { MERGE_VCF_FILES } from "./modules/local/merge_vcf_files"
include { FILTER_MERGED_VCF } from "./modules/local/filter_merged_vcf"
include { CREATE_FINAL_PLINK } from "./modules/local/create_final_plink"
include { SPLIT_INTO_CHROMOSOMES } from "./modules/local/split_into_chromosomes"
include { CREATE_REPORT } from "./modules/local/create_report"
include { GET_SAMPLES } from "./modules/local/get_samples"
include { MERGE_SAMPLES } from "./modules/local/merge_samples"

requiredParams = [
    'project',
    'output', 'chip',
    'build', 'strand_file',
    'refalt_file'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}


//load all plink files from sheet or from file pattern
if (params.input_csv != null) {
    
    plink_files = Channel.fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true, sep: ';')
        .map { row ->
            def csvDir = file(params.input_csv).parent
            def mapFile = csvDir.resolve(row.map)
            def pedFile = csvDir.resolve(row.ped)
            tuple("run_${row.run}", row.prefix, file(mapFile, checkIfExists: true), file(pedFile, checkIfExists: true))
        }

} else if (params.input != null) {

    plink_files = Channel.fromFilePairs(params.input, size: 2, flat: true, checkIfExists: true)
        .map {
            row -> tuple(row[0], row[0], row[1], row[2])
        }

} else {
    exit 1, "Parameter 'input' or 'input_csv' is required."
}

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

    vcf_files = FILTER_AND_FIX_STRAND_FLIPS.out.vcf_files
    vcf_files_index = FILTER_AND_FIX_STRAND_FLIPS.out.vcf_files_index

    if (params.reference.vcf != null) {
        study_files = vcf_files.map{
            it -> tuple(it, "study")
        }
        study_files = study_files.concat(Channel.of(tuple(file(params.reference.vcf, checkIfExists: true), "reference")))
        GET_SAMPLES (
            study_files
        )

        MERGE_SAMPLES (
            GET_SAMPLES.out.collect()
        )

        vcf_files = vcf_files.concat(Channel.of(file(params.reference.vcf, checkIfExists: true)))
        vcf_files_index = vcf_files_index.concat(Channel.of(file(params.reference.vcf  + ".tbi", checkIfExists: true)))

    }

    //TODO: write a samples.csv file: sample,type (reference or study) --> this file is then used as input to pgs-reporter

    MERGE_VCF_FILES (
        vcf_files.collect(),
        vcf_files_index.collect()
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
