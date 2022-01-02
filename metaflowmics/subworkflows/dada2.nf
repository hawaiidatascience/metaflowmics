#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.early_chimera_removal = params.early_chimera_removal ?: false

module_dir = "../modules"

include{ DADA2_FILTERANDTRIM } from "$module_dir/R/dada2/filterAndTrim/main.nf" \
    addParams( options: [publish_dir: "quality-filtering"] )
include{ DADA2_DEREPFASTQ } from "$module_dir/R/dada2/derepFastq/main.nf" \
    addParams( options: [publish_dir: "quality-filtering"] )
include{ DADA2_CHIMERA } from "$module_dir/R/dada2/removeBimeraDenovo/main.nf" \
    addParams( options: [publish_dir: "chimera"] )
include{ DADA2_LEARNERRORS } from "$module_dir/R/dada2/learnErrors/main.nf" \
    addParams( options: [publish_dir: "denoising"] )
include{ DADA2_DADA } from "$module_dir/R/dada2/dada/main.nf" \
    addParams( options: [publish_dir: "denoising"] )
include{ DADA2_MERGEPAIRS } from "$module_dir/R/dada2/mergePairs/main.nf" \
    addParams( options: [publish_dir: "read-merging"] )

// Other imports
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read-tracking"] )


workflow DADA2 {
    take:
    reads // (meta, fwd, rev) tuples

    main:
    tracked_reads = Channel.empty()

	/*
	 ========================================================================================
	 Quality filtering
	 ========================================================================================
	 */	

    qc = DADA2_FILTERANDTRIM( reads )
    derep = DADA2_DEREPFASTQ( qc.fastq ).rds

	/*
	 ========================================================================================
	 Chimera filtering (ITS pipeline design choice)
	 ========================================================================================
	 */

	// Requires single end reads or merged paired-end at this stage
    if ( params.early_chimera_removal ) {
        nochim = DADA2_CHIMERA( derep )
        tracked_reads = tracked_reads.mix(nochim.summary)
        derep = nochim.rds
    }

	/*
	 ========================================================================================
	 Error model and denoising
	 ========================================================================================
	 */

	if (params.paired_end) {
		derep = derep.multiMap {
			R1: [it[0] + [orient: "R1"], it[1][0]]
			R2: [it[0] + [orient: "R2"], it[1][1]]
		}
		derep = derep.R1.mix(derep.R2)
	} else {
		derep = derep.map{[it[0] + [orient: "R1"], it[1]]}
	}

	if ( params.pool == "F" ) {
		// error model built on all fastq files independently (N_samples * 1 or 2)
		err = DADA2_LEARNERRORS( derep )
	} else {
		// error model built on all fwd files and all reverse files (1 or 2 models)
		derep = derep.map{[[id: it[0].orient], it[1]]}.groupTuple(by: 0)
		err = DADA2_LEARNERRORS( derep )
	}

	dada = DADA2_DADA( derep.join(err.rds) )
	
	/*
	 ========================================================================================
	 Read merging
	 ========================================================================================
	 */
	
    merged = DADA2_MERGEPAIRS(
        derep.collect{it[1]},
        dada.denoised.collect{it[1]}
    )

	/*
	 ========================================================================================
	 Read tracking
	 ========================================================================================
	 */
	
    merge_summary = SUMMARIZE_TABLE(
        merged.count_table.map{[it[0], "read-merging", it[1]]}
    ).map{it[1]}

    tracked_reads = tracked_reads.mix(qc.summary).mix(dada.summary).mix(merge_summary)

    emit:
    fasta = merged.fasta
    fasta_dup = merged.fasta_dup // only for single_end 
    count_table = merged.count_table
    tracking = tracked_reads
}

