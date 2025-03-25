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
include { GATK4_MUTECT2 } from '../modules/nf-core/gatk4/mutect2/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/picard/createsequencedictionary/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFMITNANEXT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta

    main:
    // ch_samplesheet.view { i -> "samplesheet: ${i}" }
    // reference_genome.view { i -> "reference_genome: ${i}" }
    CHOPPER (ch_samplesheet,ch_fasta)
    
    
    MINIMAP2_ALIGN (
            CHOPPER.out.fastq ,
            tuple ([ id:'test_ref'], file(ch_fasta)),
            true,
            true,
            false,
            false
        )
    ch_versions = Channel.empty()

    emit:
    "Hola"


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
