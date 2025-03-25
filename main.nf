#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/nfmitnanext
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/nfmitnanext
    Website: https://nf-co.re/nfmitnanext
    Slack  : https://nfcore.slack.com/channels/nfmitnanext
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NFMITNANEXT  } from './workflows/nfmitnanext'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_nfmitnanext_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_nfmitnanext_pipeline'
//include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_nfmitnanext_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
//params.fasta = getGenomeAttribute('fasta')

params.input = ""
params.outdir = ""
params.contamination_fasta = ""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_NFMITNANEXT {

    take:
    samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_contamination_fasta
    

    main:

    // ch_samplesheet = Channel.fromPath(params.input).view(i -> 'out main first $i')


    // Channel
    // .fromPath(samplesheet)
    // .splitCsv(header: true)
    // .map { row ->
    //     def meta = [
    //         id: row.sample,
    //         single_end: true
    //     ]
    //     def fastq = file(row.fastq_1)
    //     return [meta, fastq]
    // }
    // .set { ch_samplesheet }    
    //
    // WORKFLOW: Run pipeline
    //
    NFMITNANEXT (
    //    ch_samplesheet,
        samplesheet,
        ch_fasta,
        ch_contamination_fasta
    )

    // emit:
    //     NFMITNANEXT.out // TODO: not set, but it must return something like assembly, or several data
    //     ch_samplesheet
    //     ch_fasta
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    //params.input = "/home/andresfl/NF-MITNANEX/testing_hack/samples.csv"
    //
    ch_fasta = params.fasta
    ch_contamination_fasta = params.contamination_fasta

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
    NFCORE_NFMITNANEXT (
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_fasta,
        ch_contamination_fasta
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    // PIPELINE_COMPLETION (
    //     params.email,
    //     params.email_on_fail,
    //     params.plaintext_email,
    //     params.outdir,
    //     params.monochrome_logs,
    //     params.hook_url,
    //    // NFCORE_NFMITNANEXT.out.multiqc_report
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
