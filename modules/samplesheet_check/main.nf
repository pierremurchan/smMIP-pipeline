process SAMPLESHEET_CHECK {

    input:
    path samplesheet

    output:
    path '*.csv' , emit: csv

        script:
    """
    python /home/pierre/Documents/CMD_PLACEMENT/smMIP/dev_data/bin/check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv
    """

}