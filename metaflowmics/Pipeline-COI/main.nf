#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !params.single_end

module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// module imports
include{ PEAR } from "$module_dir/pear/main.nf" \
    addParams( options: [publish_dir: "interm/read_processing/merging"] )
include{ CUTADAPT } from "$module_dir/cutadapt/main.nf" \
    addParams( options: [publish_dir: "interm/read_processing/demux"])
include{ READ_TRACKING } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )

// subworkflow imports
include { DADA2 } from "$subworkflow_dir/dada2.nf" \
    addParams( outdir: "$params.outdir/interm/read_processing",
              trunc_len: "0", trunc_quality: 2, min_read_len: 20,              
              early_chimera_removal: false, format: "mothur" )
include { MOTHUR } from "$subworkflow_dir/mothur.nf" \
	addParams( is_coding: true, outdir: "$params.outdir/interm/contig_processing" )
include { HOLOVIEWS } from "$subworkflow_dir/holoviews.nf" \
    addParams( options: [publish_dir: "figures"] )
include { DIVERSITY } from "$subworkflow_dir/diversity.nf" \
    addParams( options: [publish_dir: "postprocessing"] )

// Main workflow
workflow pipeline_COI {
    take:
    reads
    barcodes
	db

    main:

    // For read tracking through the pipeline
    read_tracking = Channel.empty()
    
	/*
	 ========================================================================================
	 Demultiplexing and read merging
	 ========================================================================================
	 */	
    demux = CUTADAPT(
        reads,
        barcodes
    )
    // Convert reads from (dataset, all_reads) to channel emitting (sample_name, reads)
    sample_reads = demux.reads.map{ it[1] }.flatten()
        .map{[it.getSimpleName().replaceFirst(/_R[12]$/, ""), it]}
        .groupTuple(by: 0)
        .map{[[id: it[0], paired_end: params.paired_end], it[1]]}

	read_tracking = read_tracking.mix(
		sample_reads.map{"demux,,${it[0].id},${it[1][0].countFastq()}"}
	)
	
    if (params.single_end || params.merge_with != "PEAR") {
        merged = sample_reads
    }
    else {
        merged = PEAR( sample_reads ).assembled.map{
            {[[id: it[0].id, paired_end: false], it[1]]}
        }
		read_tracking = read_tracking.mix(
			merged.map{"pear,,${it[0].id},${it[1].countFastq()}"}
		)
    }

	/*
	 ========================================================================================
     Read trimming, QC, denoising and merging
	 ========================================================================================
	 */	
    asvs = DADA2( merged )

    otus = MOTHUR(
        asvs.fasta,
        asvs.count_table,
		db
    )
	
	/*
	 ========================================================================================
     Read tracking through the pipeline
	 ========================================================================================
	 */

	read_tracking = read_tracking.collectFile(newLine: true)
		.mix(asvs.tracking).collectFile(name: "read_summary.csv")

    READ_TRACKING(
		otus.tracking.combine(read_tracking)
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
        otus.constaxonomy
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
    reads = Channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .map{[ [id: it[0]], it[1] ]}
    barcodes = Channel.fromPath(params.barcodes, checkIfExists: true).collect()

    // Important: DBs need to be value channels?
	dbs = Channel.fromPath("${params.db_dir}/*.{afa,tax}")
		.map{[[db_name: it.getSimpleName(),
			   db_type: it.getName().tokenize(".")[1],
			   id: it.getName()],
			  it]}

    pipeline_COI(reads, barcodes, dbs)
}
