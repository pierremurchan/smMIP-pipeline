process OUTPUT_CSV {
     publishDir "${params.outdir}/", mode: 'copy'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path (high_conf_snvs)
    path (high_conf_indels)
    path (phenotype_file)

    output:
    path 'output_summary.csv'
    path "versions.yml", emit: versions

    script:
    """
    python ${workflow.projectDir}/bin/smmip_pipeline_output_csv.py --cleaned-bams-dir ${workflow.projectDir}/${params.outdir}/smMIP-tools/cleaned_bams \\
                                                     --high-conf-snvs ${high_conf_snvs} \\
                                                     --high-conf-indels ${high_conf_indels} \\
                                                     --configuration-csv ${phenotype_file} \\
                                                     --sex-coverage-threshold 50 \\
                                                     --output output_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}