process VARIANT_REPORT {
    publishDir "${params.outdir}/smMIP-tools/", mode: 'copy'

    //conda "${moduleDir}/environment.yml"
    //container "library://murchanp/smmip-pipeline/smmiptools:sha256.f59694cb11dc219dca68f370ecf7a653f1a9348af25f54eba5c797e205425173"

    input:
    path called_mutations
    path coverage_heatmap
    path phenotype_file

    output:
    path 'variant_report.html'

    script:
    """
    python3 ${workflow.projectDir}/bin/variant_report.py --called-mutations ${called_mutations} \\
                                                         --template-path ${workflow.projectDir}/assets/ \\
                                                         --coverage-heatmap ${coverage_heatmap} \\
                                                         --configuration-csv ${phenotype_file} \\
                                                         --output-file variant_report.html

    """

    stub:
    """
    touch variant_report.html
    """
}