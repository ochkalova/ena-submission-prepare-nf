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
        .set { downloaded_reads_ch }

    downloaded_reads_ch.forward
        .groupTuple()
        .branch { meta, reads_list ->
            single: reads_list.size() == 1
            multi: reads_list.size() > 1
        }
        .set { forward_reads_ch }

    downloaded_reads_ch.reverse
        .groupTuple()
        .branch { meta, reads_list ->
            single: reads_list.size() == 1
            multi: reads_list.size() > 1
        }
        .set { reverse_reads_ch }

    CONCATENATE_FORWARD(forward_reads_ch.multi)
    CONCATENATE_REVERSE(reverse_reads_ch.multi)

    merged_forward_reads = forward_reads_ch.single.mix(CONCATENATE_FORWARD.out.file_out)
    merged_reverse_reads = reverse_reads_ch.single.mix(CONCATENATE_REVERSE.out.file_out)

    fasta
        .join(merged_forward_reads)
        .join(merged_reverse_reads)
        .multiMap { id, assembly, forward, reverse ->
            assembly: assembly
            reads: [id, [forward, reverse] ]
        }
        .set { bbmap_input_ch }

    BBMAP_ALIGN(bbmap_input_ch.reads, bbmap_input_ch.assembly)

    // emit:
    
}