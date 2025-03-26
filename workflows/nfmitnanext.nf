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
include { GATK4_FILTERMUTECTCALLS } from '../modules/nf-core/gatk4/filtermutectcalls/main'
include { BCFTOOLS_ANNOTATE } from '../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_DLOOP } from '../modules/nf-core/bcftools/annotate/main'


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
    // TODO some parameters are not set yet
    // gatk Mutect2 -R $ref_genome -L $ID --mitochondria-mode \
    // --annotation "DepthPerAlleleBySample" --min-pruning $min_pruning \
    // $kmer_size -I $aln_file -O $vcf_nofilt_file

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
    
    // ch_prueba=GATK4_MUTECT2.out.vcf.join(GATK4_MUTECT2.out.tbi).join(GATK4_MUTECT2.out.stats)
    // ch_prueba.view()
    // Filter varaints 
    //gatk FilterMutectCalls --mitochondria-mode -O $vcf_file -R $ref_genome -V $vcf_nofilt_file
    GATK4_FILTERMUTECTCALLS(
        GATK4_MUTECT2.out.vcf.join(GATK4_MUTECT2.out.tbi).join(GATK4_MUTECT2.out.stats).map{meta,vcf,tbi,stats -> tuple(meta,vcf,tbi,stats,[],[],[],[])},
        tuple([ id:'test_ref', single_end:true], file(ch_fasta)),
        tuple([ id:'test_ref', single_end:true], file(params.fai)),
        PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict,
        true
    )

    BCFTOOLS_ANNOTATE(
            GATK4_FILTERMUTECTCALLS.out.vcf.map{ meta,vcf -> tuple (meta, vcf, [], params.cds_bed, [])},
            params.header_cds,
            params.rename_chr,
            "CDS"
    )

    BCFTOOLS_ANNOTATE_DLOOP(
            BCFTOOLS_ANNOTATE.out.vcf.map{ meta,vcf -> tuple (meta, vcf, [], params.dloop_bed, [])},
            params.header_dloop,
            params.rename_chr,
            "DLOOP"
    )

    emit:
    "Hola"


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
