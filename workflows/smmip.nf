#!/usr/bin/env nextflow

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.phenotype) { ch_phenotype_file = file(params.phenotype) } else { exit 1, 'Phenotype file not specified!' }
if (params.design_file) { ch_design_file = file(params.design_file) } else { exit 1,  "No design file provided!" }
if (params.bwa) { index_ch = file(params.bwa) } else { exit 1, "No BWA index files provided!" }

ch_fasta = Channel.empty() // ch_fasta is needed by certain nf-core modules, but we will have no need for it

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
include { FASTQC_TRIMGALORE } from '../subworkflows/fastqc_trimgalore/main.nf'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/bam_sort_stats_samtools/main.nf'

include { BWA_MEM } from '../modules/bwamem/main.nf'
include { MULTIQC } from '../modules/multiqc/main.nf'
include { SMMIP_COVERAGE_HEATMAP } from '../modules/coverage_heatmap/main.nf'
include { OUTPUT_CSV } from '../modules/output_csv/main.nf'
include { VARIANT_REPORT } from '../modules/variant_report/main.nf'

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

    ch_fasta.map{ it -> [ [id: it.simpleName], [it] ]
    }.set{ fasta_tuple }

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    // Validate input samplesheet and phenotype file

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

    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    // MODULE:
    // Align reads with BWA-MEM
    
    //BWA_MEM( ch_cat_fastq, index_ch, sort_bam )
    //.bam
    //.set { ch_bam }

    BWA_MEM( FASTQC_TRIMGALORE.out.reads, index_ch, sort_bam )
    .bam
    .set { ch_bam }

    // SUBWORKFLOW:
    // Generate alignment statistics

    BAM_SORT_STATS_SAMTOOLS( BWA_MEM.out.bam, fasta_tuple )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    ch_reports = ch_reports.mix(BAM_SORT_STATS_SAMTOOLS.out.stats.map{ meta, stats -> stats})
    ch_reports = ch_reports.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.map{ meta, flagstat -> flagstat})
    ch_reports = ch_reports.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.map{ meta, idxstats -> idxstats})

    // MODULE:
    // Run MultiQC

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )

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

    // Module:
    // generate heatmap of coverage per smMIP per sample
    SMMIP_COVERAGE_HEATMAP( SMMIP_TOOLS.out.map_smmips_done.map { it[1] }.collect() )

    // Module:
    // generate output csv
    // TO DO: include variants in output csv
    OUTPUT_CSV( SMMIP_TOOLS.out.map_smmips_done.map { it[1] }.collect() )

    // Module:
    // generate a variant report
    // NOT YET IMPLEMENTED
    VARIANT_REPORT( SMMIP_TOOLS.out.called_mutations, SMMIP_COVERAGE_HEATMAP.out, ch_phenotype_file)
}