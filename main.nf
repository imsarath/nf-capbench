#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/capbench
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/capbench
    Website: https://nf-co.re/capbench
    Slack  : https://nfcore.slack.com/channels/capbench
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAPBENCH  } from './workflows/capbench'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_capbench_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_capbench_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_capbench_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.ref_genome_fasta   = getGenomeAttribute('fasta')
params.ref_genome_fai     = getGenomeAttribute('fai')
params.ref_genome_dict    = getGenomeAttribute('dict')
params.bwamem2_index      = getGenomeAttribute('bwamem2_index')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CAPBENCH {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    //
    // Initialise channels for reference genome
    //
    ch_genome_fasta  = params.ref_genome_fasta  ? Channel.fromPath(params.ref_genome_fasta).map{ it -> [[id:'genome_fasta'], it]}.collect() : Channel.empty()
    ch_genome_fai    = params.ref_genome_fai    ? Channel.fromPath(params.ref_genome_fai).map{ it -> [[id:'genome_fai'], it]}.collect() : Channel.empty()
    ch_dict          = params.ref_genome_dict   ? Channel.fromPath(params.ref_genome_dict).map{ it -> [[id:'genome_dict'], it]}.collect() : Channel.empty()
    ch_bwamem2_index = params.bwamem2_index     ? Channel.fromPath(params.bwamem2_index).map{ it -> [[id:'bwamem2_index'], it]}.collect() : Channel.empty()

    //
    // WORKFLOW: Run pipeline
    //
    CAPBENCH (
        samplesheet,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_bwamem2_index
    )
    emit:
    multiqc_report = CAPBENCH.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CAPBENCH (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_CAPBENCH.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
