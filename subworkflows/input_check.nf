//
// Check input files
//

include { SAMPLESHEET_CHECK } from '../modules/samplesheet_check/main.nf'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    //phenotype                                 // Channel: [ pheno ]
    reads                                     // channel: [ val(meta), [ reads ] ]
    //versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample

    fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    
    return fastq_meta
}

/*
def examine_phenotype(pheno){

    Channel
        .fromPath(pheno)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_cols = ['condition']

        if (!row.keySet().containsAll(expected_cols)) exit 1, "error: 'condition' is not a column name in the phenotype file.\n\nThe primary response variable must be named 'condition', please refer to the usage documentation online"

        def condition  = row.condition.matches('NA') ? 'NA' : row.condition

        if(condition == '') exit 1, "error: Invalid phenotype file, condition column contains empty cells."
        if(condition.matches('NA')) exit 1, "error: NA value in phenotype condition column."

        }
        .toList()

        return Channel.value(file(pheno))
}
*/