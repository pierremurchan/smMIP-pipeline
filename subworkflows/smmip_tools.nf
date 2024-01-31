include { MAP_SMMIPS } from '../modules/map_smmips/main.nf'
include { PILEUPS } from '../modules/pileups/main.nf'
include { CALL_MUTATIONS } from '../modules/call_mutations/main.nf'

workflow SMMIP_TOOLS {

    take:
    ch_bam_to_map
    ch_design_file
    ch_annotated_design_file
    ch_phenotype_file

    main:
    MAP_SMMIPS( ch_bam_to_map, ch_design_file )
        .clean_bam
        .set { ch_cleaned_bam }
    
    // Module:
    // Pileups
    PILEUPS( ch_cleaned_bam, ch_design_file )
        .pileups_done
        .map { it[1] }
        .collect()
        .set { ch_pileups_done }

    // Module:
    // Call mutations
    CALL_MUTATIONS ( ch_phenotype_file, ch_annotated_design_file, ch_pileups_done, ch_bam_to_map )

}