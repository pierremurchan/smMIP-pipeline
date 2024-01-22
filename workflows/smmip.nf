#!/usr/bin/env nextflow

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.phenotype) { ch_phenotype_file = file(params.phenotype) } else { exit 1, 'Phenotype file not specified!' }
if (params.bwa) { index_ch = file(params.bwa) } else { error "No BWA index files provided!" }

if (params.design_file == null) {
    error "No design file provided!"
}

//if (params.email && !params.email.matches("[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\\.[a-zA-Z]{2,6}")) {
//    error "Invalid email address provided: ${params.email}"
//}

ch_design_file = params.design_file ? file(params.design_file) : null
ch_annotated_design_file = params.annotated_design_file ? file(params.annotated_design_file) : null

design_file_channel = params.design_file ? Channel.fromPath(params.design_file) : Channel.empty()

// Set optional parameters
sort_bam = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/input_check.nf'
include { QC } from '../subworkflows/fastqc/fastqc.nf'
//include { PROCESS_READS } from '../subworkflows/process_reads.nf'
include { BWAMEM2_MEM } from '../modules/bwamem2/main.nf'
include { ANNOTATE_SNVs } from '../modules/annotate_snvs/main.nf'

//include { SMMIP_TOOLS } from '../subworkflows/smmip_tools.nf'
include { MAP_SMMIPS } from '../modules/map_smmips/main.nf'
include { PILEUPS } from '../modules/pileups/main.nf'
include { CALL_MUTATIONS } from '../modules/call_mutations/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SMMIP {

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    // Validate input samplesheet

    INPUT_CHECK( ch_input )
    .reads
    .map { meta, fastq ->
        meta.id = meta.id.split('_')[0..-2].join('_')
        [ meta, fastq ]
    }
    .set { ch_fastq }
    
    // SUBWORKFLOW: 
    // Run FastQC
    QC( ch_fastq )

    // MODULE: 
    // Align reads
    BWAMEM2_MEM(ch_fastq, index_ch, sort_bam)
    .bam
    .set { ch_bam }

    ch_bam.map { it[1] }.set { ch_bam_path }

    //
    // 2. smMIP Analysis
    //

    // Module:
    // Annotate design file
    if (params.annotated_design_file && file(params.annotated_design_file).exists()) {
        annotated_design_file = Channel.fromPath(params.annotated_design_file)
    } else {
        annotated_design_file = ANNOTATE_SNVs( ch_design_file )
    }
    
    // Module:
    // Map smMIPs
    // will integrate into subworkflow
    MAP_SMMIPS( ch_bam, ch_design_file )
    .set { ch_cleaned_bam }

    // Module:
    // Pileups
    // will integrate into subworkflow
    PILEUPS( ch_cleaned_bam, ch_design_file )
    .pileups_done
    .map { it[1] }
    .set { ch_pileups_done }

    ch_pileups_done.view()

    // Module:
    // Call mutations
    // will integrate into subworkflow
    CALL_MUTATIONS ( ch_phenotype_file, ch_annotated_design_file, ch_pileups_done, ch_bam )
}