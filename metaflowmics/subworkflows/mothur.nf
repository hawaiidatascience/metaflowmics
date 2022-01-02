#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
module_dir = "../modules"
interm_dir = "interm/contig_processing"

// workflows used by mothur main workflow
include { MSA } from "./msa.nf" \
     addParams( outdir: "$params.outdir/$interm_dir/msa-filter" )
include { TAXONOMY } from "./taxonomy.nf" \
     addParams( outdir: "$params.outdir/$interm_dir/taxonomy" )
include { CONSENSUS as CONSENSUS_PRELIM } from "./mothur-util" \
    addParams( outdir: "$params.outdir/raw" )    
include { CONSENSUS } from "./mothur-util" \
     addParams( outdir: "$params.outdir/results" )
include { SYNC } from "./mothur-util" \
     addParams( outdir: "$params.outdir/results" )

// modules
include { MOTHUR_SCREEN_SEQS } from "$module_dir/mothur/screenSeqs/main.nf" \
    addParams( options: [publish_dir: "$interm_dir/msa-filter"] )
include { MOTHUR_CHIMERA } from "$module_dir/mothur/chimera/main.nf" \
    addParams( options: [publish_dir: "$interm_dir/chimera-filter"] )
include { MOTHUR_CLUSTER } from "$module_dir/mothur/cluster/main.nf" \
    addParams( options: [publish_dir: "$interm_dir/clustering"] )
include { MOTHUR_REMOVE_LINEAGE } from "$module_dir/mothur/removeLineage/main.nf" \
    addParams( options: [publish_dir: "$interm_dir/lineage-filter"] )
include { MOTHUR_REMOVE_RARE } from "$module_dir/mothur/removeRare/main.nf" \
    addParams( options: [publish_dir: "$interm_dir/rare-otu-filter"] )
include { GET_SUBSAMPLING_THRESHOLD } from "$module_dir/util/misc/main.nf"
include { MOTHUR_SUBSAMPLE } from "$module_dir/mothur/subsample/main.nf" \
    addParams( options: [publish_dir: "$interm_dir/subsampling"] )
include { MOTHUR_DIST_SEQS } from "$module_dir/mothur/distSeqs/main.nf" \
    addParams( cutoff: 1-params.lulu_min_match/100, format: "vsearch" )
include { LULU } from "$module_dir/R/lulu/main.nf" \
    addParams( options: [publish_dir: "$interm_dir/lulu-filter"] )
include { MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: "nseqs-sobs" )
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf"


workflow MOTHUR {
    take:
    fasta
    count_table
	dbs

    main:
	/*
	 ========================================================================================
	 Multiple sequence alignment: subworkflow also used for ITS and COI
	 ========================================================================================
	 */
	msa = MSA( 	// align.seqs()
		fasta,
		count_table,
		dbs
	)

    msa_filt = MOTHUR_SCREEN_SEQS(
        msa.aln.combine(msa.count_table)
    )

	/*
	 ========================================================================================
	 Chimeric contig filtering
	 ========================================================================================
	 */	
    chimera = MOTHUR_CHIMERA(
        msa_filt.fasta.join(msa_filt.count_table)
    )
	/*
	 ========================================================================================
	 Taxonomic assignment: subworkflow used for ITS and COI
	 ========================================================================================
	 */	
    tax = TAXONOMY( // classify.seqs()
        chimera.fasta,
		chimera.count_table,
		dbs
    ).taxonomy

	/*
	 ========================================================================================
	 OTU clustering with all the provided identity threshold
	 ========================================================================================
	 */
    cluster = MOTHUR_CLUSTER(
        chimera.fasta.join(chimera.count_table).join(tax),
        params.clustering_thresholds.split(",").collect{it as int}
    )

	/*
	 ========================================================================================
	 Prepare variables for optional tasks. We overwrite them each time we have an optional
	 input executed.
	 + add OTU identity to all channels
	 ========================================================================================
	 */
    fasta = cluster.fasta
    count_table = cluster.count_table
    taxonomy = cluster.taxonomy
	list = cluster.list
	shared = cluster.shared

	/*
	 ========================================================================================
	 Intermediate results before all optional steps to provide a first glance at the data
	 ========================================================================================
	 */
    consensus_raw = CONSENSUS_PRELIM(
        fasta,
        count_table,
        taxonomy,
        list,
        shared
    )

	// prepare read tracking channel for optional steps
	// tracked = shared.map{[[step: "clustering", otu_id: it[0].otu_id], it[1]]}
	tracked = shared.map{[it[0], "clustering", it[1]]}

	/*
	 ========================================================================================
	 Optional: remove contigs matching specific taxonomic keywords 
	 (e.g. mitochondria, chloroplasts) or not matching (e.g. unknown kingdom)
	 ========================================================================================
	 */
    if ( params.taxa_to_filter != "" ) {
        (count_table, taxonomy, list, shared) = MOTHUR_REMOVE_LINEAGE(
            count_table.join(taxonomy).join(list)
        )
        // tracked = tracked.mix(shared.map{[[step: "taxa-filter", otu_id: it[0].otu_id], it[1]]})
		tracked = tracked.mix(shared.map{[it[0], "taxa-filter", it[1]]})		
    }

	/*
	 ========================================================================================
	 Optional: discard rare contigs
	 ========================================================================================
	 */
    (list, shared, count_table) = MOTHUR_REMOVE_RARE(
        list.join(count_table)
    )
    // tracked = tracked.mix(shared.map{[[step: "rare-otus-filter", otu_id: it[0].otu_id], it[1]]})
	tracked = tracked.mix(shared.map{[it[0], "rare-otus-filter", it[1]]})		

	/*
	 ========================================================================================
	 Optional: subsampling
	 ========================================================================================
	 */
    if (!params.skip_subsampling) {
        subsampling_level = GET_SUBSAMPLING_THRESHOLD( count_table.first().map{it[1]} )
        (shared, list, count_table) = MOTHUR_SUBSAMPLE(
            list.join(count_table),
            subsampling_level
        )
        // tracked = tracked.mix(shared.map{[[step: "subsampling", otu_id: it[0].otu_id], it[1]]})
        tracked = tracked.mix(shared.map{[it[0], "subsampling", it[1]]})
    }

	
	/*
	 ========================================================================================
	 Optional: co-occurrence pattern correction with LULU
	 This step is very slow and should be skipped for large datasets (> ~1000 samples)
	 ========================================================================================
	 */
    if (!params.skip_lulu) {
        dists = MOTHUR_DIST_SEQS( consensus_raw.repfasta ).dist
        (shared, repfasta, discard, summary) = LULU(
            dists.join(shared).join(consensus_raw.repfasta)
        )
        tracked = tracked.mix(shared.map{[it[0], "lulu", it[1]]})
    }

	/*
	 ========================================================================================
	 Subset OTUs to make them the same between files
	 (since filtering is not always done on all the files each time)
	 ========================================================================================
	 */	
    results = SYNC(
        fasta.mix(count_table).mix(taxonomy),
        list,
        shared
    )

	/*
	 ========================================================================================
	 Consensus to provide:
	 - representative OTU sequence 
	 - consensus taxonomy
	 - database file if the option is set
	 ========================================================================================
	 */	
    consensus = CONSENSUS(
        results.fasta,
        results.count_table,
        results.taxonomy,
        results.list,
        shared
    )

	/*
	 ========================================================================================
	 Count tracking through the mothur subworkflow
	 ========================================================================================
	 */	
    // count tables
	tracking_ct = SUMMARIZE_TABLE(
		msa_filt.count_table.map{[it[0], "MSA", it[1]]}.mix(
			chimera.count_table.map{[it[0], "chimera-filter", it[1]]}
		)
	)

    // shared
	tracking_shared = MOTHUR_SUMMARY_SINGLE(
		tracked
	).summary

	tracking = tracking_ct.mix(tracking_shared)
		.collectFile(){["${it[0].db_name}_${it[0].db_type}.csv", it[1]]}
		.map{[it.getSimpleName(), it]}

    emit:
    repfasta = consensus.repfasta
    constaxonomy = consensus.constaxonomy
    shared = shared
    tracking = tracking
}
