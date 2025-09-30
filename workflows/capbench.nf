/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_capbench_pipeline'
include { SAMTOOLS_INDEX as INDEX_RAWBAM        } from '../modules/nf-core/samtools/index/main'


include { DOWNSAMPLING_FACTOR                                } from '../modules/local/downsampling_factor/main'
include { PICARD_BEDTOINTERVALLIST                           } from '../modules/nf-core/picard/bedtointervallist/main'
include { ALIGNMENT                                          } from '../subworkflows/local/fastq_align_bwamem2/main.nf'
include { BAM_QC_METRICS as BAMQC                            } from '../subworkflows/local/bam_qc_metrics/main'

include { BAM_QC_METRICS as BAMQC_FOR_UMIRAWBAM              } from '../subworkflows/local/bam_qc_metrics/main'
include { BAM_QC_METRICS as BAMQC_FOR_UMIPROCESSEDBAM        } from '../subworkflows/local/bam_qc_metrics/main'
include { BAM_QC_METRICS as BAMQC_FOR_DOWNSAMPLING           } from '../subworkflows/local/bam_qc_metrics/main'
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO as UMI_PROCESSING } from '../subworkflows/nf-core/fastq_create_umi_consensus_fgbio/main'
include { DOWNSAMPLING                                       } from '../subworkflows/local/downsampling/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CAPBENCH {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_genome_fasta
    ch_genome_fai
    ch_dict
    ch_bwamem2_index
    ch_panel_bed
    ch_interval_list

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run FastP
    //
    FASTP (
        ch_samplesheet,
        [], // adapter_fasta: not used in this pipeline
        params.discard_trimmed_pass,
        params.save_trimmed_fail,
        params.save_merged
    )

    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_input_reads = FASTP.out.reads

    if (!params.interval_list && params.panel_bed) {

        PICARD_BEDTOINTERVALLIST (
            ch_panel_bed,
            ch_dict,
            []
        )

        ch_versions = ch_versions.mix(PICARD_BEDTOINTERVALLIST.out.versions.first())
        ch_interval_list = PICARD_BEDTOINTERVALLIST.out.intervallist
    }

    if (params.umi_structure) {

        ch_input_reads
            .map { meta, reads ->
                return [meta.id, [meta , reads]]
            }
            .groupTuple()
            .map { id, grouped_reads ->
                def metas = grouped_reads.collect{it[0]}
                def files = grouped_reads.collect{it[1]}.flatten()
                return [metas[0], files]
            }
            .set { ch_input_reads }

        CAT_FASTQ (
            ch_input_reads
        )

        UMI_PROCESSING(
            CAT_FASTQ.out.reads,
            [],
            ch_genome_fasta,
            ch_bwamem2_index,
            ch_dict,
            "paired",
            "bwa-mem2",
            params.duplex,
            params.min_reads,
            params.min_baseq,
            params.max_base_error_rate,
            ch_interval_list.collect { it[1] },
            true
        )

        INDEX_RAWBAM (
            UMI_PROCESSING.out.mappedbam
        )

        SAMTOOLS_INDEX(
            UMI_PROCESSING.out.mappedconsensusbam
        )

        ch_umi_raw_bam        = UMI_PROCESSING.out.mappedbam
            .join(INDEX_RAWBAM.out.bai)
        ch_umi_processed_bam  = UMI_PROCESSING.out.mappedconsensusbam
            .join(SAMTOOLS_INDEX.out.bai)

        ch_versions = ch_versions.mix(UMI_PROCESSING.out.versions.first())

        BAMQC_FOR_UMIRAWBAM (
            ch_umi_raw_bam,
            ch_genome_fasta,
            ch_genome_fai,
            ch_dict,
            ch_interval_list,
            false
        )

        BAMQC_FOR_UMIPROCESSEDBAM (
            ch_umi_processed_bam,
            ch_genome_fasta,
            ch_genome_fai,
            ch_dict,
            ch_interval_list,
            false
        )

        ch_multiqc_files = ch_multiqc_files.mix(BAMQC_FOR_UMIRAWBAM.out.multiple_metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAMQC_FOR_UMIRAWBAM.out.coverage_metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAMQC_FOR_UMIPROCESSEDBAM.out.multiple_metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAMQC_FOR_UMIPROCESSEDBAM.out.coverage_metrics.collect{it[1]}.ifEmpty([]))

        ch_versions = ch_versions.mix(BAMQC_FOR_UMIRAWBAM.out.versions.first())
        ch_versions = ch_versions.mix(BAMQC_FOR_UMIPROCESSEDBAM.out.versions.first())

        ch_duplex_metrics = UMI_PROCESSING.out.duplex_metrics

        ch_duplex_metrics
            .map { meta, duplex_metrics -> duplex_metrics }
            .collect()
            .map { files ->
                def meta = [id: 'duplex_metrics']
                return([meta, files])
            }.set{ ch_duplex_metrics }


        DOWNSAMPLING (
            ch_duplex_metrics,
            ch_umi_raw_bam, // input BAM files from alignment workflow
            params.downsampling_targets,
            params.downsampling_seed,
            ch_genome_fasta,
            ch_bwamem2_index,
            ch_dict,
            ch_interval_list.collect { it[1] },
        )


        BAMQC_FOR_DOWNSAMPLING (
            DOWNSAMPLING.out.bam,
            ch_genome_fasta,
            ch_genome_fai,
            ch_dict,
            ch_interval_list,
            false
        )

    } else {

        ALIGNMENT(
            ch_input_reads,
            ch_genome_fasta,
            ch_genome_fai,
            ch_bwamem2_index
        )

        ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.dedup_metrics.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions.first())
        ch_aligned_bam = ALIGNMENT.out.dedup_bam
            .join(ALIGNMENT.out.dedup_bai)

        BAMQC (
            ch_aligned_bam,
            ch_genome_fasta,
            ch_genome_fai,
            ch_dict,
            ch_interval_list,
            true
        )

        ch_multiqc_files = ch_multiqc_files.mix(BAMQC.out.multiple_metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAMQC.out.coverage_metrics.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(BAMQC.out.versions.first())

    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'capbench_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
