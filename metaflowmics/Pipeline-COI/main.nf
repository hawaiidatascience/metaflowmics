#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !(params.single_end)
    
module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// Modules
include{ DOWNLOAD_IBOL } from "$module_dir/util/download/main.nf" \
    addParams( db_release: params.ibol_release, level: params.translation_table_mapping_level )
// include{ GET_TRANSLATION_TABLES } from "$module_dir/util/entrez-direct/main.nf" \
//     addParams( level: params.translation_table_mapping_level )
// include{ TRANSLATE } from "$module_dir/util/biopython/translation/main.nf" \
//     addParams( level: params.translation_table_mapping_level )
// include{ READ_TRACKING } from "$module_dir/util/misc/main.nf" \
//     addParams( options: [publish_dir: "read_tracking"] )

// Subworkflows
include { msa_agg as msa_per_species } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/species"], sep: ',', field: 7 )
include { msa_agg as msa_per_genus } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/genus"], sep: ',', field: 6 )
include { msa_agg as msa_per_subfamily } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/subfamily"], sep: ',', field: 5 )
include { msa_agg as msa_per_family } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/family"], sep: ',', field: 4 )
include { msa_agg as msa_per_order } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/order"], sep: ',', field: 3 )
include { msa_agg as msa_per_class } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/class"], sep: ',', field: 2 )
include { msa_agg as msa_per_phylum } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/phylum"], sep: ',', field: 1 )
include { msa_agg as msa_per_kingdom } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa/kingdom"], sep: ',', field: 0 )

include { msa_update as msa_update_per_species } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa_updates/species"], sep: ',', field: 0 )
include { msa_update as msa_update_per_genus } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa_updates/genus"], sep: ',', field: 6 )
include { msa_update as msa_update_per_subfamily } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa_updates/subfamily"], sep: ',', field: 5 )
include { msa_update as msa_update_per_family } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa_updates/family"], sep: ',', field: 4 )
include { msa_update as msa_update_per_order } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa_updates/order"], sep: ',', field: 3 )
include { msa_update as msa_update_per_class } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa_updates/class"], sep: ',', field: 2 )
include { msa_update as msa_update_per_phylum } from "$subworkflow_dir/msa.nf" \
    addParams( options: [publish_dir: "msa_updates/phylum"], sep: ',', field: 1 )

// Functions
include { helpMessage ; saveParams } from "./util.nf"


// Main workflow
workflow pipeline_COI {
    take:
    fasta

    main:
    // Download iBOL db
    versions = Channel.fromList(
        (2..25).collect{String.format('%.2f', it/4)}
    )

    if (!params.test) {
        db = DOWNLOAD_IBOL(versions)
        db_tsv = db.tsv.collectFile(keepHeader: true, skip: 1, storeDir: params.outdir)
    } else {
        db = Channel.fromPath("../../tests/COI/iBOL_COI_sub.tsv")
            .splitCsv(header: true, sep: '\t')
            .map{">$it.seqentryid;k__Animalia,p__$it.phylum_reg,c__$it.class_reg,o__$it.order_reg,f__$it.family_reg,sf__$it.subfamily_reg,g__$it.genus_reg,s__$it.species_reg\n$it.aminoraw"}
            .collectFile(newLine: true)
            .branch{faa: true}
        db_tsv = Channel.fromPath("../../tests/COI/iBOL_COI_sub.tsv")
    }
    
    // MSA per taxonomic group
    aln_species = msa_per_species(db.faa)
    aln_genus = msa_per_genus(aln_species.repr)
    aln_subfamily = msa_per_subfamily(aln_genus.repr)
    aln_family = msa_per_family(aln_subfamily.repr)
    aln_order = msa_per_order(aln_family.repr)
    aln_class = msa_per_class(aln_order.repr)
    aln_phylum = msa_per_phylum(aln_class.repr)
    aln_kingdom = msa_per_kingdom(aln_phylum.repr)

    // Realign the original sequences against the representative sequences
    aln_phylum = msa_update_per_phylum(1, aln_phylum.afa, aln_kingdom.afa)
    aln_class = msa_update_per_class(2, aln_class.afa, aln_phylum.afa)
    aln_order = msa_update_per_order(3, aln_order.afa, aln_class.afa)
    aln_family = msa_update_per_family(4, aln_family.afa, aln_order.afa)
    aln_subfamily = msa_update_per_subfamily(5, aln_subfamily.afa, aln_family.afa)
    aln_genus = msa_update_per_genus(6, aln_genus.afa, aln_subfamily.afa)    
    aln_species = msa_update_per_species(7, aln_species.afa, aln_genus.afa)

    // Translate contigs to proteins
    // translation_tables = GET_TRANSLATION_TABLES(db.taxa)
    //     .collectFile(newLine: true)    

    // asv_aa = TRANSLATE(fasta, translation_tables)
}

workflow {
    // reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1)
    //     .map{[ [id: it[0]], it[1] ]}
    // saveParams()
    pipeline_COI(Channel.fromPath(params.fasta))
}
