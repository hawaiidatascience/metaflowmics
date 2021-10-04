#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.is_coding = false

module_dir = "../modules"

include { MOTHUR_CLASSIFY_SEQS } from "$module_dir/mothur/classifySeqs/main.nf"
include{ DADA2_ASSIGN_TAXONOMY } from "$module_dir/R/dada2/assignTaxonomy/main.nf"
// include { VSEARCH_SINTAX } from "$module_dir/vsearch/taxonomy/main.nf"
include { RDP_CLASSIFY } from "$module_dir/rdp/classify/main.nf"


def map_combine(channel_a, channel_b, key){
    channel_a
        .map{ it -> [it[0][key], it] }
        .combine(channel_b.map{it -> [it[0][key], it[1..-1]]}, by: 0)
        .map { it[1] + it[2] }
}

workflow TAXONOMY {
    take:
    contigs
	count_table
    dbs

	main:
	/*
	 ========================================================================================
	 split databases
	 ========================================================================================
	 */
	dbs = dbs
		.filter{it[0].db_type != "prot"}
		.map{[[db_name: it[0].db_name], it[1]]}
		.branch{
		dada2: (it[0].db_name =~/unite/) // for ITS pipeline
		rdp: (it[0].db_name =~ /rdp/) // for COI pipeline?
		mothur: true // for 16S and COI pipeline
	}

	/*
	 ========================================================================================
     RDP with mothur (with alignment)
	 ========================================================================================
	 */

	taxonomy_mothur = MOTHUR_CLASSIFY_SEQS(
		map_combine(contigs.join(count_table), dbs.mothur.groupTuple(), "db_name")
	).taxonomy

	/*
	 ========================================================================================
     RDP with dada2(without alignment)
	 ========================================================================================
	 */
	taxonomy_dada2 = DADA2_ASSIGN_TAXONOMY(
		contigs,
		dbs.dada2
	).taxonomy

	/*
	 ========================================================================================
     RDP (without alignment)
	 ========================================================================================
	 */
	taxonomy_rdp = Channel.empty()
	// TODO

	emit:
	taxonomy = taxonomy_mothur.mix(taxonomy_dada2).mix(taxonomy_rdp)
}
