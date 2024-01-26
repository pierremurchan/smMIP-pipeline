process MAP_SMMIPS {
    publishDir "${params.outdir}/smMIP-tools/cleaned_bams/${meta.id}", 
               mode: 'copy'

    input:
    tuple val(meta), path(bamFile)
    path design_file

    output:
    tuple val(meta), path("${meta.id}_clean.bam"), emit: clean_bam
    tuple val(meta), path("${meta.id}_clean.bam.bai")
    tuple val(meta), path("${meta.id}_clean.sam")
    tuple val(meta), path("${meta.id}_filtered.sam")
    tuple val(meta), path("${meta.id}_filtered_read_counts.txt")
    tuple val(meta), path("${meta.id}_UMI_usage_per_smMIP.txt")
    tuple val(meta), path("${meta.id}_raw_coverage_per_smMIP.txt")

    script:
    """
    Rscript ${workflow.projectDir}/bin/map_smMIPs_extract_UMIs.R --bam.file ${bamFile} \
                                                    --panel.file ${design_file} \
                                                    --sample.name ${meta.id} \
                                                    --output \${PWD} \
                                                    --filtered.reads ${params.filtered_reads} \
                                                    --OVERLAP ${params.OVERLAP} \
                                                    --MAPQ ${params.MAPQ} \
                                                    --threads ${params.threads} \
                                                    --code ${workflow.projectDir}/bin/
    """

    stub:
    """
    touch  ${meta.id}_clean.bam
    touch  ${meta.id}_clean.bam.bai
    touch  ${meta.id}_clean.sam
    touch  ${meta.id}_filtered_read_counts.txt
    touch  ${meta.id}_raw_coverage_per_smMIP.txt
    touch  ${meta.id}_UMI_useage_per_smMIP.txt
    """
}
