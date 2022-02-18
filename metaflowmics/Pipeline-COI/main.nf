#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !params.single_end

module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// module imports
include{ PEAR } from "$module_dir/pear/main.nf" \
    addParams( options: [publish_dir: "interm/read_processing/merging"] )
include{ CUTADAPT; CUTADAPT_JAMP } from "$module_dir/cutadapt/main.nf" \
    addParams( options: [publish_dir: "interm/read_processing/demux"])
include{ RDP_CLASSIFY } from "$module_dir/rdp/classify/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/taxonomy"] )
include{ MOTHUR_CHIMERA } from "$module_dir/mothur/chimera/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/chimera-filter"] )
include { MOTHUR_CLUSTER } from "$module_dir/mothur/cluster/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/clustering"] )
include { MOTHUR_REMOVE_LINEAGE } from "$module_dir/mothur/removeLineage/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/lineage-filter"] )
include { MOTHUR_REMOVE_RARE } from "$module_dir/mothur/removeRare/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/rare-otu-filter"] )
include { GET_SUBSAMPLING_THRESHOLD } from "$module_dir/util/misc/main.nf"
include { MOTHUR_SUBSAMPLE } from "$module_dir/mothur/subsample/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/subsampling"] )
include { MOTHUR_CONSENSUS as MOTHUR_CONSENSUS_RAW } from "$module_dir/mothur/consensus/main.nf" \
    addParams( options: [publish_dir: "raw/consensus"] )    
include { MOTHUR_CONSENSUS } from "$module_dir/mothur/consensus/main.nf" \
     addParams( options: [publish_dir: "results/consensus"] )
include { MOTHUR_SYNC } from "$module_dir/mothur/sync/main.nf" \
    addParams( options: [publish_dir: "results/detail"] )
include { MOTHUR_MAKE_DATABASE as MOTHUR_MAKE_DATABASE_RAW } from "$module_dir/mothur/makeDatabase/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_MAKE_DATABASE } from "$module_dir/mothur/makeDatabase/main.nf" \
    addParams( options: [publish_dir: "results"] )

include { MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: "nseqs-sobs" )
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf"
include{ READ_TRACKING } from "$module_dir/util/misc/main.nf" \
	addParams( options: [publish_dir: "read_tracking"] )

// subworkflow imports
include { DADA2 } from "$subworkflow_dir/dada2.nf" \
    addParams( outdir: "$params.outdir/interm/read_processing",
              // trunc_len: "0,0", trunc_quality: 2, min_read_len: 20, paired_end: false,
              early_chimera_removal: false, format: "mothur" )
include { HOLOVIEWS } from "$subworkflow_dir/holoviews.nf" \
    addParams( options: [publish_dir: "figures"] )
include { DIVERSITY } from "$subworkflow_dir/diversity.nf" \
    addParams( options: [publish_dir: "postprocessing"], skip_unifrac: true )

// Main workflow
workflow pipeline_COI {
    take:
    reads
    barcodes
	db

    main:

    // For read tracking through the pipeline
    tracking_reads = Channel.empty()
    
	/*
	 ========================================================================================
	 Demultiplexing and read merging
	 ========================================================================================
	 */

	if (!params.skip_demux) {
		demux = params.jamp_demux ?
			CUTADAPT_JAMP(reads, barcodes) :
			CUTADAPT(reads, barcodes)
		
		// Convert reads from (dataset, all_reads) to channel emitting (sample_name, reads)
		sample_reads = demux.reads.map{ it[1] }.flatten()
			.map{[it.getSimpleName().replaceFirst(/_R[12]$/, ""), it]}
			.groupTuple(by: 0)
			.map{[[id: it[0], paired_end: params.paired_end], it[1]]}
	} else {
		sample_reads = reads
	}
		
	tracking_reads = tracking_reads.mix(
		sample_reads.map{"demux,,${it[0].id},${it[1][0].countFastq()}"}
	).collectFile(newLine: true)
	
    // merged = PEAR( sample_reads ).assembled
	// tracking_reads = tracking_reads.mix(
	// 	merged.map{"pear,,${it[0].id},${it[1].countFastq()}"}
	// ).collectFile(newLine: true)
	
	/*
	 ========================================================================================
     Read trimming, QC, denoising and merging
	 ========================================================================================
	 */

	asvs = DADA2( sample_reads )
    // asvs = DADA2( merged )
	chimera = MOTHUR_CHIMERA(
        asvs.fasta.join( asvs.count_table )
    )
	tax = RDP_CLASSIFY( chimera.fasta, db ).taxonomy
	
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
    consensus_raw = MOTHUR_CONSENSUS_RAW(
		shared.join(list).join(fasta).join(count_table).join(taxonomy)
    )

	if (params.compute_mothur_db) {
		MOTHUR_MAKE_DATABASE_RAW(
			shared
				.join(consensus_raw.constaxonomy)
				.join(consensus_raw.repfasta)
				.join(consensus_raw.repcount_table)
		)
	}

	// prepare read tracking channel for optional steps
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
	tracked = tracked.mix(shared.map{[it[0], "rare-otus-filter", it[1]]})		

	/*
	 ========================================================================================
	 Optional: subsampling
	 ========================================================================================
	 */

	if (!params.skip_subsampling) {
		if (params.custom_subsampling_level > 0) {
			subsampling_level = Channel.value(params.custom_subsampling_level)
		} else {
			subsampling_level = GET_SUBSAMPLING_THRESHOLD( count_table.first().map{it[1]} )
		}
		(shared, list, count_table) = MOTHUR_SUBSAMPLE(
			list.join(count_table),
			subsampling_level
		)
        tracked = tracked.mix(shared.map{[it[0], "subsampling", it[1]]})
    }

	/*
	 ========================================================================================
	 Subset OTUs to make them the same between files
	 (since filtering is not always done on all the files each time)
	 ========================================================================================
	 */	

	sync = MOTHUR_SYNC(
		shared.join(list).join(fasta).join(count_table).join(taxonomy)
	)

	/*
	 ========================================================================================
	 Consensus to provide:
	 - representative OTU sequence 
	 - consensus taxonomy
	 - database file if the option is set
	 ========================================================================================
	 */	

	consensus = MOTHUR_CONSENSUS(
		shared.join(sync.list).join(sync.fasta).join(sync.count_table).join(sync.taxonomy)
	)

	if (params.compute_mothur_db) {
		MOTHUR_MAKE_DATABASE(
			shared
				.join(consensus.constaxonomy)
				.join(consensus.repfasta)
				.join(consensus.repcount_table)
		)
	}

	/*
	 ========================================================================================
	 Read tracking through the pipeline
	 ========================================================================================
	 */

	// reads
	tracking_reads = tracking_reads
		.mix(asvs.tracking)
		.collectFile(name: "tracking.csv")

	// count tables
	tracking_ct = SUMMARIZE_TABLE(
		chimera.count_table.map{[it[0], "chimera-filter", it[1]]}
	).map{it[1]}

    // shared
	tracking_shared = MOTHUR_SUMMARY_SINGLE(
		tracked
	).summary.map{it[1]}
	
	tracking = READ_TRACKING(
		tracking_reads.mix(tracking_ct).mix(tracking_shared)
			.collectFile(name: "tracking_all.csv")
	)
	
	/*
	 ========================================================================================
	 Interactive visualization with holoviews python package
	 ========================================================================================
	 */	
    HOLOVIEWS(
        shared,
        consensus.constaxonomy,
    )

	/*
	 ========================================================================================
	 Diversity metrics (alpha, beta + phylogenetic tree)
	 ========================================================================================
	 */	
    DIVERSITY(
        consensus.repfasta,
        shared
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .map{[ [id: it[0], paired_end: params.paired_end], it[1] ]}
    barcodes = Channel.fromPath(params.barcodes).collect()

	db = Channel.fromPath("${params.db_dir}/*.{txt,xml,properties}").collect()

    pipeline_COI(reads, barcodes, db)
}
