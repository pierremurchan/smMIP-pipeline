#!/usr/bin/env nextflow

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.phenotype) { ch_phenotype_file = file(params.phenotype) } else { exit 1, 'Phenotype file not specified!' }

bwa_index = Channel.fromPath(params.bwa)
if (params.bwa) { index_ch = file(params.bwa) } else { error "No BWA index files provided!" }
//if (params.fasta) { fasta_ch = file(params.fasta)} else { error "No FASTA file provided!" }


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
include { SMMIP_TOOLS } from '../subworkflows/smmip_tools.nf'

include { CAT_FASTQ } from '../modules/cat_fastq/main.nf'
include { FASTQC } from '../modules/fastqc/main.nf'
//include { BWAMEM2_MEM } from '../modules/bwamem2/main.nf'
include { BWA_MEM } from '../modules/bwamem/main.nf'
include { MULTIQC } from '../modules/multiqc/main.nf'
//include { SAMTOOLS_INDEX } from '../modules/samtools_index/main.nf'
//include { PICARD_COLLECTHSMETRICS } from '../modules/collecthsmetrics/main.nf'
include { SMMIP_COVERAGE_HEATMAP } from '../modules/coverage_heatmap/main.nf'

include { ANNOTATE_SNVs } from '../modules/annotate_snvs/main.nf'
include { MAP_SMMIPS } from '../modules/map_smmips/main.nf'
include { PILEUPS } from '../modules/pileups/main.nf'
include { CALL_MUTATIONS } from '../modules/call_mutations/main.nf'

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

    INPUT_CHECK( ch_input )
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
        annotated_design_file = Channel.fromPath(params.annotated_design_file)
    } else {
        annotated_design_file = ANNOTATE_SNVs( ch_design_file )
    }

    ch_bam_to_map = ch_input_files.is_bam.mix( ch_bam )

    // Subworkflow:
    // smMIP-tools process
    SMMIP_TOOLS( ch_bam_to_map,  ch_design_file, ch_annotated_design_file, ch_phenotype_file )

    SMMIP_COVERAGE_HEATMAP( SMMIP_TOOLS.out.map_smmips_done.map { it[1] }.collect() )
}