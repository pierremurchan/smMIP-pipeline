include { MAP_SMMIPS } from '../modules/map_smMIPs/main.nf'
//include { PILEUPS } from '../modules/pileups/main.nf'
//include { CALL_MUTATIONS } from '../modules/call_mutations/main.nf'

workflow SMMIP_TOOLS {
    take:
    path bamFile
    path annotated_design_file

    main:
    bamFiles_ch = Channel.fromPath("${params.bam_dir}/*.bam")

}