process MAP_SMMIPS {
    publishDir "${params.outdir}/smMIP-tools/cleaned_bams", 
               mode: 'copy'

    //conda "${moduleDir}/environment.yml"
    container "library://murchanp/smmip-pipeline/smmiptools:sha256.f59694cb11dc219dca68f370ecf7a653f1a9348af25f54eba5c797e205425173"

    input:
    tuple val(meta), path(bamFile)
    path design_file

    output:
    tuple val(meta), path("${meta.id}_clean.bam"), emit: clean_bam
    tuple val(meta), path("${meta.id}_clean.bam.bai")
    //tuple val(meta), path("${meta.id}_clean.sam")
    tuple val(meta), path("${meta.id}_filtered_read_counts.txt")
    tuple val(meta), path("${meta.id}_raw_coverage_per_smMIP.txt")
    tuple val(meta), path("${meta.id}_UMI_usage_per_smMIP.txt")
    tuple val(meta), path("tmp/map_smmips_${meta.id}.done.tmp"), emit: map_smmips_done

    script:
    """
    Rscript ${workflow.projectDir}/bin/map_smMIPs_extract_UMIs.R --bam.file ${bamFile} \
                                                    --panel.file ${design_file} \
                                                    --sample.name ${meta.id} \
                                                    --output \${PWD} \
                                                    --filtered.reads y \
                                                    --OVERLAP ${params.OVERLAP} \
                                                    --MAPQ ${params.MAPQ} \
                                                    --threads ${params.threads} \
                                                    --code ${workflow.projectDir}/bin/
    mkdir -p tmp
    touch tmp/map_smmips_${meta.id}.done.tmp
    """

    stub:
    """
    mkdir -p tmp
    touch  ${meta.id}_clean.bam
    touch  ${meta.id}_clean.bam.bai
    #touch  ${meta.id}_clean.sam
    touch  ${meta.id}_filtered_read_counts.txt
    touch  ${meta.id}_raw_coverage_per_smMIP.txt
    touch  ${meta.id}_UMI_usage_per_smMIP.txt
    touch tmp/map_smmips_${meta.id}.done.tmp
    """
}
