#!/usr/bin/env nextflow

include { FASTQC } from '../../modules/fastqc/main.nf'
include { MULTIQC } from '../../modules/multiqc/main.nf'

workflow QC {
    take:
    reads

    main:
    
    FASTQC( reads ).set { fastqc_ch }
    MULTIQC( fastqc_ch.collect() )
    
    }