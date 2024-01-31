process QC_PLOTS {
    publishDir "${params.outdir}/smMIP-tools/", mode: 'copy'

    input:
    path pileupsDone

    output:
    path "Cohort_Coverage_per_smMIP.txt"

    script:
    """
    Rscript ${workflow.projectDir}/bin/QC_plots.R --dir ${workflow.projectDir}/${params.outdir}/smMIP-tools/pileups \
                                                  --code ${workflow.projectDir}/bin 
    """

    stub:
    """
    touch Cohort_Coverage_per_smMIP.txt
    """
}