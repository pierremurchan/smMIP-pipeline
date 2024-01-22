process ANNOTATE_SNVs {
    publishDir "${params.outdir}/smMIP-tools", mode:'copy'

    input:
    path design_file
    //path outdir

    output:
    path "annotated_${design_file.baseName}.txt"


    script:
    """
    Rscript ${workflow.projectDir}/bin/Annotate_SNVs.R --panel.file ${design_file} \
                                          --threads ${params.threads} \
                                          --genome ${params.genome} \
                                          --species "hsapiens" \
                                          --threads ${params.threads} \
                                          --code ${workflow.projectDir}/bin/
    """

    stub:
    """
    touch "annotated_${design_file.baseName}.txt"
    """
}