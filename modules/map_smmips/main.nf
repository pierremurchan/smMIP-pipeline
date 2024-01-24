process MAP_SMMIPS {
    publishDir "${params.outdir}/smMIP-tools/cleaned_bams/", 
               mode: 'copy'

    input:
    tuple val(meta), path(bamFile)
    path design_file

    output:
    tuple val(meta), path("${meta.id}/*_clean.bam"), emit: clean_bam

    script:
    """
    #mkdir -p ${meta.id}
    Rscript ${workflow.projectDir}/bin/map_smMIPs_extract_UMIs.R --bam.file ${bamFile} \
                                                    --panel.file ${design_file} \
                                                    --sample.name ${meta.id} \
                                                    --output \${PWD} \
                                                    --filtered.reads ${params.filtered_reads} \
                                                    --OVERLAP ${params.OVERLAP} \
                                                    --MAPQ ${params.MAPQ} \
                                                    --threads ${params.threads} \
                                                    --code ${workflow.projectDir}/bin/
    #cp ${meta.id}/* ./
    #mv ${meta.id}/* ./
    #rm -r ${meta.id}
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
