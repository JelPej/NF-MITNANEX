/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//include { FASTQC                 } from '../modules/nf-core/fastqc/main'
//include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nfmitnanext_pipeline'
include { CHOPPER } from '../modules/nf-core/chopper/main'
include { MINIMAP2_ALIGN } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_VIEW } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_BAM2FQ } from '../modules/nf-core/samtools/bam2fq/main'
include { FLYE } from '../modules/nf-core/flye/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFMITNANEXT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_contamination_fasta
    ch_min_mapQ
    ch_flye_mode

    main:
    // ch_samplesheet.view { i -> "samplesheet: ${i}" }
    // reference_genome.view { i -> "reference_genome: ${i}" }
    CHOPPER (ch_samplesheet,ch_contamination_fasta)
    
    
    MINIMAP2_ALIGN (
            CHOPPER.out.fastq ,
            tuple ([ id:'test_ref'], file(ch_fasta)),
            true,
            "bai",
            false,
            false
        )

    ch_versions = Channel.empty()

    // def index = MINIMAP2_ALIGN.out.index.last()
    // index.view()
    // ch_samtools_input = MINIMAP2_ALIGN
    //     .out
    //     .bam
    //     .view()
    //     .map { tuple ->
    //             def meta = tuple[0]
    //             def input = tuple[1]
    //             //def index = MINIMAP2_ALIGN.out.index.last()
    //             return tuple(meta, input,index)
    //         }
    //     //.view()
    // ch_samtools_input.view()

    ch_samtools_input =  MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index).view()
    SAMTOOLS_VIEW (
        ch_samtools_input,
        tuple ([ id:'test_ref'], file(ch_fasta)),
        ch_min_mapQ)
    
    SAMTOOLS_BAM2FQ(
        SAMTOOLS_VIEW.out.bam,
        false
    )
    //ch_flye_mode.view()
    FLYE (
        SAMTOOLS_BAM2FQ.out.reads,
        ch_flye_mode
    )

    //FLYE(ch_samplesheet,ch_contamination_fasta)
    // FLYE_ALIGN (
    //     SAMTOOLS.out.bam , 
    //     tuple([])
    // )


    emit:
    "Hola"


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
