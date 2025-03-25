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
<<<<<<< HEAD
//params.fasta = getGenomeAttribute('fasta')
=======
params.fasta = getGenomeAttribute('fasta')
//params.samplesheet = "/home/andresfl/NF-MITNANEX/testing_hack/samples.csv"
params.input = "/home/andresfl/NF-MITNANEX/testing_hack/samples.csv"
params.outdir = 'results'

// Channel.fromPath(params.samplesheet).view(i -> 'out main param def $i')
>>>>>>> master

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
<<<<<<< HEAD
        samplesheet // channel: samplesheet read in from --input
        reference_genome // string: Path to reference genome 
=======
    samplesheet // channel: samplesheet read in from --input
    ch_fasta
    
>>>>>>> master

    main:

    // ch_samplesheet = Channel.fromPath(params.input).view(i -> 'out main first $i')


    Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row ->
        def meta = [
            id: row.sample,
            single_end: true
        ]
        def fastq = file(row.fastq_1)
        return [meta, fastq]
    }
    .set { ch_samplesheet }    
    ch_fasta = params.fasta
    //
    // WORKFLOW: Run pipeline
    //
    NFMITNANEXT (
<<<<<<< HEAD
        samplesheet,
        reference_genome 
    )

    emit:
    NFMITNANEXT.out // TODO: not set, but it must return something like assembly, or several data
=======
        ch_samplesheet,
        ch_fasta
    )
//     emit:
//     multiqc_report = NFMITNANEXT.out.multiqc_report // channel: /path/to/multiqc_report.html
>>>>>>> master
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
<<<<<<< HEAD
        params.reference_genome
=======
        params.genome
>>>>>>> master
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
<<<<<<< HEAD
    //     params.hook_url,
    //    // NFCORE_NFMITNANEXT.out.multiqc_report
=======
    //     params.hook_url
    //     // NFCORE_NFMITNANEXT.out.multiqc_report
>>>>>>> master
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
