process CALL_MUTATIONS {
    publishDir "${params.outdir}/smMIP-tools/", mode: 'copy'

    //conda "${moduleDir}/environment.yml"
    container "library://murchanp/smmip-pipeline/smmiptools:sha256.f59694cb11dc219dca68f370ecf7a653f1a9348af25f54eba5c797e205425173"


    input:
    path phenotype_file
    path anno_design_file
    path pileupsDone
    tuple val(meta), path(bamFile)

    output:
    path "called_mutations.txt", emit: mutation_calls

    script:
    """
    Rscript ${workflow.projectDir}/bin/calling_mutations.R --summary ${workflow.projectDir}/${params.outdir}/smMIP-tools/pileups \
                                                           --file ${phenotype_file} \
                                                           --alleles ${anno_design_file} \
                                                           --binomial ${params.binomial} \
                                                           --overlap.coverage ${params.overlap_coverage} \
                                                           --maf ${params.maf} \
                                                           --vaf ${params.vaf} \
                                                           --pval ${params.pval} \
                                                           --drop ${params.drop} \
                                                           --output \${PWD} \
                                                           --threads ${params.threads} \
                                                           --code ${workflow.projectDir}/bin 
    """

    stub:
    """
    touch called_mutations.txt
    """
}