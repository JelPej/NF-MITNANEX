process MITOMASTER_SUBMIT {
    tag "$meta.id"
    label 'process_low'
    
    input:
        tuple val (meta), path(vcf)
    output:
        path("*_mitomaster_result.txt")
    script:
    """
    curl -s -X POST https://mitomap.org/mitomaster/websrvc.cgi \\
        -F "file=@${vcf}" \\
        -F "fileType=snvlist" \\
        -F "output=detail" \\
        -o ${meta.id}_mitomaster_result.txt
    """
}
