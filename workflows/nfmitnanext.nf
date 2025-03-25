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
    
    // Map reads to reference genome
    MINIMAP2_ALIGN (
            CHOPPER.out.fastq ,
            tuple ([ id:'test_ref'], file(ch_fasta)),
            true,
            "bai",
            false,
            false
        )
    
    Channel
        .fromPath(ch_fasta)
        .map { fasta_file ->
        def meta = [ id: fasta_file.getBaseName().replaceAll(/\.(fa|fasta)(\.gz)?$/, '') ]
        tuple(meta, fasta_file)
        }
        .set { ch_dict_name }
        
    // Create dict for variant calling
    PICARD_CREATESEQUENCEDICTIONARY (
        ch_dict_name
    )

    // Variant calling for mitochondria

    GATK4_MUTECT2(
        MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index).map{meta,bam,bai -> tuple(meta,bam,bai,[])},
        tuple ([ id:'genome'], file(ch_fasta)),
        tuple([ id:'genome'],file(params.fai)),
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict,
        [],
        [],
        [],
        [],
        true
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
