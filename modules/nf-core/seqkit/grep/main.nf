process SEQKIT_GREP {
    tag "$meta.id"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(sequence)
    //path pattern
    tuple val(meta2), path(flye_asm_info)

    output:
    tuple val(meta), path("*.{fa,fq}.gz")  , emit: filter
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def suffix = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"
    //def pattern_file = pattern ? "-f ${pattern}" : ""
    def flye_asm_info = "${flye_asm_info}"
    """
    mito_contig_ID=\$(awk -F'\t' '\$4 ~ /Y/ { if (\$3 > max_val) { max_val = \$3; result = \$1 } } END { print result }' ${flye_asm_info})
    echo \$mito_contig_ID
    seqkit \\
        grep \\
        $args \\
        --threads $task.cpus \\
        -p \$mito_contig_ID \\
        ${sequence} \\
        -o ${prefix}.${suffix}.gz \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def suffix = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"

    """
    echo "" | gzip > ${prefix}.${suffix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
