include { SRATOOLS_PREFETCH } from "../../modules/nf-core/sratools/prefetch/main"
include { SRATOOLS_FASTERQDUMP } from "../../modules/nf-core/sratools/fasterqdump/main"
include { CAT_CAT as CONCATENATE_FORWARD } from "../../modules/nf-core/cat/cat/main"
include { CAT_CAT as CONCATENATE_REVERSE } from "../../modules/nf-core/cat/cat/main"
include { BBMAP_ALIGN } from '../../modules/nf-core/bbmap/align/main'

workflow  COMPUTE_ASSEMBLY_COVERAGE {
    take:
    fasta      // [ meta, path(fasta) ]
    reads      // [ meta, tuple(reads_acc) ]

    main:
    SRATOOLS_PREFETCH(
        reads.transpose(),
        [],   // ncbi_settings
        [],   // certificate
    )

    SRATOOLS_FASTERQDUMP(
        SRATOOLS_PREFETCH.out.sra,
        [],   // ncbi_settings
        [],   // certificate
    )

    SRATOOLS_FASTERQDUMP.out.reads
        .multiMap { id, reads_accessions ->
            forward: [id, reads_accessions[0] ]
            reverse: [id, reads_accessions[1] ]
        }
        .set { cat_input_ch }
    
    CONCATENATE_FORWARD(cat_input_ch.forward.groupTuple())
    CONCATENATE_REVERSE(cat_input_ch.reverse.groupTuple())

    fasta
        .join(CONCATENATE_FORWARD.out.file_out)
        .join(CONCATENATE_REVERSE.out.file_out)
        .multiMap { id, assembly, forward_concatenated, reverse_concatenated ->
            assembly: [id, assembly ]
            reads: [id, [forward_concatenated, reverse_concatenated] ]
        }
        .set { bbmap_input_ch }

    BBMAP_ALIGN(bbmap_input_ch.assembly, bbmap_input_ch.reads)

    // emit:
    
}