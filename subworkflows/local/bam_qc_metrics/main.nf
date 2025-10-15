
include { PICARD_COLLECTMULTIPLEMETRICS  } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS        } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_COLLECTWGSMETRICS       } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { JUMBLE_RUN                     } from '../../../modules/local/jumble/run/main'

workflow BAM_QC_METRICS {
    take:
    ch_input_bam // channel: input BAM files from alignment workflow
    ch_genome_fasta
    ch_genome_fai
    ch_genome_dict
    ch_interval_list
    ch_jumble_ref  // channel: reference for jumble, only required if running
    low_pass_wgs

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Collect BAM metrics with Picard
    //

    PICARD_COLLECTMULTIPLEMETRICS(
        ch_input_bam,
        ch_genome_fasta,
        ch_genome_fai
    )

    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    ch_coverage_metrics = Channel.empty()

    if (low_pass_wgs) {

        PICARD_COLLECTWGSMETRICS(
            ch_input_bam,
            ch_genome_fasta,
            ch_genome_fai,
            []
        )

        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first())
        coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTWGSMETRICS.out.metrics)

    } else {
        ch_bam = ch_input_bam
            .combine(ch_interval_list)
            .map { meta, bam, bai, meta2, interval_list ->
                [meta, bam, bai, interval_list, interval_list]
            }

        //
        // MODULE: Collect HS metrics with Picard
        //

        PICARD_COLLECTHSMETRICS(
            ch_bam,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dict,
            [[], []]
        )

        coverage_metrics = ch_coverage_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics)
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    }

    JUMBLE_RUN(
        ch_input_bam,
        ch_jumble_ref
    )


    emit:
    multiple_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    coverage_metrics = coverage_metrics
    versions         = ch_versions

}
