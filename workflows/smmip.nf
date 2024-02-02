#!/usr/bin/env nextflow

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.phenotype) { ch_phenotype_file = file(params.phenotype) } else { exit 1, 'Phenotype file not specified!' }

bwa_index = Channel.fromPath(params.bwa)
if (params.bwa) { index_ch = file(params.bwa) } else { error "No BWA index files provided!" }

if (params.design_file == null) {
    error "No design file provided!"
}

ch_design_file = params.design_file ? file(params.design_file) : null

// Set optional parameters
sort_bam = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/input_check.nf'
include { SMMIP_TOOLS } from '../subworkflows/smmip_tools.nf'

include { CAT_FASTQ } from '../modules/cat_fastq/main.nf'
include { FASTQC } from '../modules/fastqc/main.nf'
include { BWA_MEM } from '../modules/bwamem/main.nf'
include { MULTIQC } from '../modules/multiqc/main.nf'
include { SMMIP_COVERAGE_HEATMAP } from '../modules/coverage_heatmap/main.nf'

include { ANNOTATE_SNVs } from '../modules/annotate_snvs/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def multiqc_report = []

workflow SMMIP {

    ch_reports  = Channel.empty()
    ch_versions = Channel.empty()

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    // Validate input samplesheet

    INPUT_CHECK( ch_input, ch_phenotype_file )
    .reads
    .groupTuple(by: [0])
    .branch {
        meta, reads ->
            is_bam  : meta.is_bam == true
                return [ meta, reads.flatten() ]
            single_run_fastq: meta.is_bam == false && reads.size() == 1
                return [ meta, reads.flatten() ]
            multiple_runs_fastq: meta.is_bam == false && reads.size() > 1
                return [ meta, reads.flatten() ]
    }
    .set { ch_input_files }

    // MODULE:
    // Concatenate FastQ files from same sample if required
    CAT_FASTQ (
        ch_input_files.multiple_runs_fastq
    )
    .reads
    .mix(ch_input_files.single_run_fastq )
    .set{ ch_cat_fastq }

    // MODULE: 
    // Run FastQC on FASTQ files

    if (!params.skip_fastqc) {
        FASTQC( ch_cat_fastq )
        ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }
    
    BWA_MEM( ch_cat_fastq, index_ch, sort_bam )
    .bam
    .set { ch_bam }

    //
    // 2. smMIP Analysis
    //

    // Module:
    // Annotate design file
    if (params.annotated_design_file && file(params.annotated_design_file).exists()) {
        ch_annotated_design_file = Channel.fromPath(params.annotated_design_file)
    } else {
        ch_annotated_design_file = ANNOTATE_SNVs( ch_design_file )
    }

    ch_bam_to_map = ch_input_files.is_bam.mix( ch_bam )

    // Subworkflow:
    // smMIP-tools process
    SMMIP_TOOLS( ch_bam_to_map,  ch_design_file, ch_annotated_design_file, ch_phenotype_file )

    // Modeule:
    // generate heatmap of coverage per smMIP per sample
    SMMIP_COVERAGE_HEATMAP( SMMIP_TOOLS.out.map_smmips_done.map { it[1] }.collect() )
}