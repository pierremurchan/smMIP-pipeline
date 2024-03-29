process CATEGORISE_MUTATIONS {
    publishDir "${params.outdir}/smMIP-tools/", mode: 'copy'

    //conda "${moduleDir}/environment.yml"
    container "library://murchanp/smmip-pipeline/smmiptools:sha256.f59694cb11dc219dca68f370ecf7a653f1a9348af25f54eba5c797e205425173"


    input:
    path(mutation_calls)

    output:
    path "High_Confidence_Indels.txt"
    path "Lower_Confidence_Indels.txt"
    path "High_Confidence_SNVs.txt"
    path "Lower_Confidence_SNVs.txt"

    script:
    """
    Rscript ${workflow.projectDir}/bin/mutation_categories.R --input ${mutation_calls} \
                                                             --code ${workflow.projectDir}/bin
    """

    stub:
    """
    touch "High_Confidence_Indels.txt"
    touch "Lower_Confidence_Indels.txt"
    touch "High_Confidence_SNVs.txt"
    touch "Lower_Confidence_SNVs.txt"
    """
}