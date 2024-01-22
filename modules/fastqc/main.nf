process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc", mode:'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path("fastqc_${meta.id}_logs")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Reads: ${reads}"
    mkdir -p fastqc_${meta.id}_logs
    fastqc -o fastqc_${meta.id}_logs -f fastq -q ${reads.join(' ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "fastqc_${meta.id}_logs"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
