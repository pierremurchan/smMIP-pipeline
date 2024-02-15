process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode:'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path('*')

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """

    stub:
    """
    touch multiqc_report.html
    """
}
