process MERGE_SAMPLES {
    
    publishDir "$params.output", mode: 'copy'

    input:
    path(sample_files)

    output:
    path("samples.txt")

    script:
    """
    echo "sample;source" > samples.txt
    cat ${sample_files} >> samples.txt
    """

}