include { samplesheetToList    } from 'plugin/nf-schema'

include { COMPUTE_ASSEMBLY_COVERAGE } from '../subworkflows/compute_assembly_coverage/main'

workflow PIPELINE {

    println params.input

    samplesheet_ch = Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    

    input_ch = samplesheet_ch
        .multiMap { id, metagenome, runs ->
            assembly: [ id, metagenome ]
            reads: [ id, tuple(runs.split(";")) ]
        }

    COMPUTE_ASSEMBLY_COVERAGE(input_ch.assembly, input_ch.reads)

}