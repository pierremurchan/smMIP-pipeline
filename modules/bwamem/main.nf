process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/bwa", mode:'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:9c0128851101dafef65cef649826d2dbe6bedd7e-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:9c0128851101dafef65cef649826d2dbe6bedd7e-0' }"

    input:
    tuple val(meta), path(reads)
    path(index)
    val   sort_bam

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    //path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools $samtools_command $args2 --threads $task.cpus -o ${prefix}.bam -

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}