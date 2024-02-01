process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode:'copy'

    input:
    path('*')

    output:
    path('multiqc_report.html')

    script:
    """
    /miniconda/bin/multiqc .
    """

    stub:
    """
    touch multiqc_report.html
    """
}
