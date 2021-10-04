#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !(params.single_end)
    
module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// DADA2 module imports
include { DADA2 } from "$subworkflow_dir/dada2.nf" \
    addParams( outdir: "$params.outdir/interm/read_processing",
              early_chimera_removal: false, format: "mothur" )
// mothur module imports
include { MOTHUR } from "$subworkflow_dir/mothur.nf"

// Other imports
include{ DOWNLOAD_SILVA_FOR_MOTHUR } from "$module_dir/bash/download/main.nf" \
    addParams( db_release: params.silva_db )
include{ READ_TRACKING } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )

// Subworkflows
include { HOLOVIEWS } from "$subworkflow_dir/holoviews.nf" \
    addParams( options: [publish_dir: "figures"] )
include { DIVERSITY } from "$subworkflow_dir/diversity.nf" \
    addParams( options: [publish_dir: "postprocessing"] )

// Functions
include { helpMessage ; saveParams } from "./util.nf"


// Main workflow
workflow pipeline_16S {
    take:
    reads

    main:
	/*
	 ========================================================================================
     Read trimming, QC, denoising and merging
	 ========================================================================================
	 */	
	asvs = DADA2( reads )

	/*
	 ========================================================================================
	 Contigs curation with mothur (subworkflow)
	 ========================================================================================
	 */	
	
    // Download SILVA db for mothur
    db = DOWNLOAD_SILVA_FOR_MOTHUR()

    otus = MOTHUR(
        asvs.fasta,
        asvs.count_table,
        db.tax.mix(db.align)
    )
	
	/*
	 ========================================================================================
     Read tracking through the pipeline
	 ========================================================================================
	 */
	
    READ_TRACKING(
		otus.tracking.combine(
			reads.map{"raw,,${it[0].id},${it[1][0].countFastq()}"}.collectFile(newLine: true)
				.mix(asvs.tracking).collectFile(name: "summary.csv")
		)
			.map{[it[0], it[1..-1]]}.transpose()
			.collectFile(){[it[0], it[1]]}
			.map{[it.getSimpleName(), it]}
    )

	
	/*
	 ========================================================================================
	 Interactive visualization with holoviews python package
	 ========================================================================================
	 */	
    HOLOVIEWS(
        otus.shared,
        otus.constaxonomy,
    )

	/*
	 ========================================================================================
	 Diversity metrics (alpha, beta + phylogenetic tree)
	 ========================================================================================
	 */	
    DIVERSITY(
        otus.repfasta,
        otus.shared
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1)
        .map{[ [id: it[0], paired_end: params.paired_end], it[1] ]}
    saveParams()
    pipeline_16S(reads)
}
