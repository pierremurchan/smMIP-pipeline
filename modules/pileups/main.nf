process PILEUPS {
    publishDir "${params.outdir}/smMIP-tools/pileups/", mode: 'copy'

    input:
    tuple val(meta), path(cleaned_bamFile)
    path design_file_path

    output:
    tuple val(meta), path("tmp/pileups_${meta.id}.done.tmp"), emit: pileups_done // outputting tmp file to ensure that Pileups is completed before CallingMutations runs, tried using state dependency but it didn't work
    //path "${params.outdir}/smMIP-tools/pileups/", emit: pileups_dir
    tuple val(meta), path("*_raw_pileup.txt"), emit: raw_pileup
    tuple val(meta), path("*_sscs_pileup.txt"), emit: sscs_pilep

    script:
    """
    Rscript ${workflow.projectDir}/bin/smMIP_level_raw_and_consensus_pileups.R --bam.file ${cleaned_bamFile} \
                                                                               --panel.file ${design_file_path} \
                                                                               --sample.name ${meta.id} \
                                                                               --output \${PWD} \
                                                                               --mnd ${params.mnd} \
                                                                               --mmq ${params.mmq} \
                                                                               --mbq ${params.mbq} \
                                                                               --rank ${params.rank} \
                                                                               --family.size ${params.family_size} \
                                                                               --consensus.cutoff ${params.consensus_cutoff} \
                                                                               --threads ${params.threads} \
                                                                               --code ${workflow.projectDir}/bin
    mkdir -p tmp
    touch tmp/pileups_${meta.id}.done.tmp
    """

    stub:
    """
    mkdir -p tmp

    touch tmp/pileups_${meta.id}.done.tmp

    touch \${PWD}/${meta.id}_raw_pileup.txt
    touch \${PWD}/${meta.id}_sscs_pileup.txt
    """
}