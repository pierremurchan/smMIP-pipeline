process OUTPUT_CSV {
     publishDir "${params.outdir}/", mode: 'copy'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path (map_smmips_done)

    output:
    path 'coverage.csv'
    path "versions.yml", emit: versions

    script:
    """
    python ${workflow.projectDir}/bin/calculate_coverage.py --input-dir ${workflow.projectDir}/${params.outdir}/smMIP-tools/cleaned_bams \\
                                                     --sex-coverage-threshold 50 \\
                                                     --output coverage.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}