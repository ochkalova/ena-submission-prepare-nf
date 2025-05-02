process MODULE {
    label 'process_single'
    container ''

    input:
    tuple val(meta), path(assembly_fasta)

    output:
    path('file.txt'), emit: some_output

    script:
    """

    """
}