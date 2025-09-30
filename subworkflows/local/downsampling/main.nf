

include { DOWNSAMPLING_FACTOR  } from '../../../modules/local/downsampling_factor/main'
include { SAMTOOLS_SUBSAMPLE   } from '../../../modules/local/subsample/main'
include { SAMTOOLS_INDEX       } from '../../../modules/nf-core/samtools/index/main'

include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO as UMI_PROCESSING } from '../../../subworkflows/nf-core/fastq_create_umi_consensus_fgbio/main'


workflow DOWNSAMPLING {
    take:
    ch_duplex_metrics
    ch_input_bam // channel: input BAM files from alignment workflow
    targets_str
    seed
    ch_genome_fasta
    ch_bwamem2_index
    ch_dict
    interval_list

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Calculate downsampling factors
    //

    DOWNSAMPLING_FACTOR (
        ch_duplex_metrics,
        targets_str
    )


    //
    // MODULE: Downsample BAM files with samtools view -s
    //

    ch_input_subsample = ch_input_bam
        .combine(DOWNSAMPLING_FACTOR.out.factors_tsv)
        .flatMap { meta, bam, bai, factors ->
            def lines = factors.readLines().drop(1)
            lines.collect { line ->
                def fields = line.split("\t")
                def targetM = fields[3].toFloat()
                def fraction = fields[4].toFloat()
                def newmeta = meta.clone() as Map
                newmeta.id = "${meta.id}_${targetM}M"

                tuple(newmeta, bam, bai, fraction, targetM)
            }
        }
        .unique { it[0].id }

    SAMTOOLS_SUBSAMPLE (
        ch_input_subsample,
        seed
    )

    ch_downsampled_bams = SAMTOOLS_SUBSAMPLE.out.bam

    UMI_PROCESSING(
        [],
        ch_downsampled_bams,
        ch_genome_fasta,
        ch_bwamem2_index,
        ch_dict,
        "paired",
        "bwa-mem2",
        params.duplex,
        params.min_reads,
        params.min_baseq,
        params.max_base_error_rate,
        interval_list,
        false
    )

    SAMTOOLS_INDEX(
        UMI_PROCESSING.out.mappedconsensusbam
    )

    ch_downsampled_bam = UMI_PROCESSING.out.mappedconsensusbam
            .join(SAMTOOLS_INDEX.out.bai)
    //
    // COLLECT VERSIONS
    //

    ch_versions = ch_versions.mix(DOWNSAMPLING_FACTOR.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_SUBSAMPLE.out.versions.first())
    ch_versions = ch_versions.mix(UMI_PROCESSING.out.versions.first())

    emit:
    bam            = ch_downsampled_bam
    duplex_metrics = UMI_PROCESSING.out.duplex_metrics
    ch_versions
}
