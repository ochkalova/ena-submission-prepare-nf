include { SRATOOLS_PREFETCH } from "../../modules/nf-core/sratools/prefetch/main"
include { SRATOOLS_FASTERQDUMP } from "../../modules/nf-core/sratools/fasterqdump/main"
include { CAT_CAT as CONCATENATE_FORWARD } from "../../modules/nf-core/cat/cat/main"
include { CAT_CAT as CONCATENATE_REVERSE } from "../../modules/nf-core/cat/cat/main"
include { BBMAP_ALIGN } from '../../modules/nf-core/bbmap/align/main'

workflow  COMPUTE_ASSEMBLY_COVERAGE {
    take:
    fasta                       // [ meta, path(fasta) ]
    assembly_reads_mapping      // [ meta, tuple(reads_acc) ]

    main:
    
    assembly_reads_mapping
        .transpose()
        .map { meta, reads_acc ->
            [ [id:reads_acc], reads_acc ]
        }
        .unique()
        .set { for_download_ch }
    
    SRATOOLS_PREFETCH(
        for_download_ch,
        [],   // ncbi_settings
        [],   // certificate
    )

    SRATOOLS_FASTERQDUMP(
        SRATOOLS_PREFETCH.out.sra,
        [],   // ncbi_settings
        [],   // certificate
    )

    assembly_reads_mapping                            // [ [[id:ASSEMBLY_1], [SRR1, SRR2]] ]
        .transpose()                                  // [ [[id:ASSEMBLY_1], SRR1], [[id:ASSEMBLY_1], SRR2] ]
        .map { meta, read_acc ->
            [ [id:read_acc], meta ]                   // [ [[id:SRR1], [id:ASSEMBLY_1]], [[id:SRR2], [id:ASSEMBLY_1]] ]
        }
        .groupTuple()
        .mix(SRATOOLS_FASTERQDUMP.out.reads)          // [ [[id:SRR1], [id:ASSEMBLY_1]], [[id:SRR2], [id:ASSEMBLY_1]], [[id:SRR1], [SRR1_1.fastq.gz, SRR1_2.fastq.gz]], [[id:SRR2], [SRR2_1.fastq.gz, SRR2_2.fastq.gz]] ]
        .groupTuple()                                 // [ [[id:SRR1], [[id:ASSEMBLY_1], [SRR1_1.fastq.gz, SRR1_2.fastq.gz]]], [[id:SRR2], [[id:ASSEMBLY_1], [SRR2_1.fastq.gz, SRR2_2.fastq.gz]]] ]
        .map { meta, data -> 
            def assembly_list = data[0]
            def reads_files = data[1]
            [assembly_list, reads_files] }
        .transpose(by:0)
        .multiMap { meta, reads ->
            forward: [ meta, reads[0] ]               // [ [[id:ASSEMBLY_1], SRR1_1.fastq.gz], [[id:ASSEMBLY_1], SRR2_1.fastq.gz] ]
            reverse: [ meta, reads[1] ]               // [ [[id:ASSEMBLY_1], SRR1_2.fastq.gz], [[id:ASSEMBLY_1], SRR2_2.fastq.gz] ]
        }
        .set { downloaded_reads_ch }

    // This code below was used when I mainly aligned reads to NOT co-assemblies. 
    // SRATOOLS_PREFETCH(
    //     assembly_reads_mapping.transpose(),
    //     [],   // ncbi_settings
    //     [],   // certificate
    // )

    // SRATOOLS_FASTERQDUMP(
    //     SRATOOLS_PREFETCH.out.sra,
    //     [],   // ncbi_settings
    //     [],   // certificate
    // )

    // SRATOOLS_FASTERQDUMP.out.reads
    //     .multiMap { meta, reads ->
    //         forward: [meta, reads[0] ]
    //         reverse: [meta, reads[1] ]
    //     }
    //     .set { downloaded_reads_ch }

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

    merged_forward_reads = forward_reads_ch.single.transpose().mix(CONCATENATE_FORWARD.out.file_out)
    merged_reverse_reads = reverse_reads_ch.single.transpose().mix(CONCATENATE_REVERSE.out.file_out)

    fasta
        .join(merged_forward_reads)
        .join(merged_reverse_reads)
        .multiMap { meta, assembly, forward, reverse ->
            assembly: assembly
            reads: [meta, [forward, reverse] ]
        }
        .set { bbmap_input_ch }

    BBMAP_ALIGN(bbmap_input_ch.reads, bbmap_input_ch.assembly)

    // emit:
    
}