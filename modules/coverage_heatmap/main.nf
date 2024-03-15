process SMMIP_COVERAGE_HEATMAP {
    publishDir "${params.outdir}/smMIP-tools/", mode: 'copy'
    
    container "library://murchanp/smmip-pipeline/python3.8_plotly_pandas.sif:sha256.2835d23990d1960b1fa6d85b27b502c967ae3e71373945c4c8e8b3b8523b82ec"

    input:
    path (map_smmips_done)

    output:
    path 'coverage_heatmap.html'

    script:
    """
    python3 ${workflow.projectDir}/bin/coverage_heatmap.py \\
        --input-dir ${workflow.projectDir}/${params.outdir}/smMIP-tools/cleaned_bams \\
        --output coverage_heatmap.html
    """

    stub:
    """
    touch coverage_heatmap.html
    """
}