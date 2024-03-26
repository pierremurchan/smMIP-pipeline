include { MAP_SMMIPS } from '../modules/map_smmips/main.nf'
include { PILEUPS } from '../modules/pileups/main.nf'
include { CALL_MUTATIONS } from '../modules/call_mutations/main.nf'
include { CATEGORISE_MUTATIONS } from '../modules/categorise_mutations/main.nf'

workflow SMMIP_TOOLS {

    take:
    ch_bam_to_map
    ch_design_file
    ch_annotated_design_file
    ch_phenotype_file

    main:
    // Module:
    // map smMIPS
    MAP_SMMIPS( ch_bam_to_map, ch_design_file )
    ch_cleaned_bam = MAP_SMMIPS.out.clean_bam
    
    // Module:
    // generate pileups
    PILEUPS( ch_cleaned_bam, ch_design_file )
        .pileups_done
        .map { it[1] }
        .collect()
        .set { ch_pileups_done }

    // Module:
    // call mutations
    CALL_MUTATIONS ( ch_phenotype_file, ch_annotated_design_file, ch_pileups_done, ch_bam_to_map )

    // Module:
    // categorise mutations into high and low calls
    CATEGORISE_MUTATIONS( CALL_MUTATIONS.out.called_mutations )
    
    emit:
    map_smmips_done = MAP_SMMIPS.out.map_smmips_done // again, probably a better way to implement state dependecy
    called_mutations = CALL_MUTATIONS.out.called_mutations
    high_conf_indels = CATEGORISE_MUTATIONS.out.high_conf_indels
    high_conf_snvs = CATEGORISE_MUTATIONS.out.high_conf_snvs
}