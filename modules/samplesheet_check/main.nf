process SAMPLESHEET_CHECK {

    input:
    path samplesheet

    output:
    path '*.csv' , emit: csv

        script:
    """
    python ${workflow.projectDir}/bin/check_samplesheet_v2.py \\
        $samplesheet \\
        samplesheet.valid.csv
    """

}