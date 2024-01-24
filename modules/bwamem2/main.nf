process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/bwamem2", mode:'copy'

    input:
    tuple val(meta), path(reads)
    //tuple val(meta2), path(index)
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
    INDEX=\$(find -L . -name "*.amb" | sed 's/.amb//')

    /miniconda/bin/bwa-mem2 \\
        mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | /miniconda/bin/samtools $samtools_command $args2 -@ $task.cpus -o ${prefix}.bam -

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
