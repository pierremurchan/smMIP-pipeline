//
// Check input files
//

include { SAMPLESHEET_CHECK } from '../modules/samplesheet_check/main.nf'

workflow INPUT_CHECK {
    take:
    samplesheet
    ch_phenotype

    main:
    phenotype = params.phenotype ? examine_phenotype(ch_phenotype) : Channel.empty()

    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_input_file_channel(it) }
        .set { reads }

    emit:
    //phenotype                                 // Channel: [ pheno ]
    reads                                     // channel: [ val(meta), [ reads ] ]
    //versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ] or [meta [ bam ]]
def create_input_file_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    // Check the file extension to determine if it's BAM or FASTQ
    meta.is_bam = row.bam ? true : false

    if (meta.is_bam) {
        input_file_meta = [meta, [file(row.bam)]]
    } else {
        input_file_meta = [meta, [file(row.fastq_1), file(row.fastq_2)]]
    }
    
    return input_file_meta
}


def examine_phenotype(pheno){

    Channel
        .fromPath(pheno)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_cols = ['id', 'type', 'replicate']

        if (!row.keySet().containsAll(expected_cols)) exit 1, "error: unexpected column name format in the phenotype file, columns should be 'id', 'type', and 'replicate'."

        def type  = row.type.matches('NA') ? 'NA' : row.type

        if(type == '') exit 1, "error: Invalid phenotype file, type column contains empty cells."
        if(type.matches('NA')) exit 1, "error: NA value in phenotype type column."

        }
        .toList()

        return Channel.value(file(pheno))
}
