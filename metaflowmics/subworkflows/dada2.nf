#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.early_chimera_removal = params.early_chimera_removal ?: false

module_dir = "../modules"

include{ DADA2_FILTERANDTRIM } from "$module_dir/R/dada2/filterAndTrim/main.nf" \
    addParams( options: [publish_dir: "quality-filtering"] )
include{ DADA2_DEREPFASTQ } from "$module_dir/R/dada2/derepFastq/main.nf" \
    addParams( options: [publish_dir: "chimera"] )
include{ DADA2_LEARNERRORS } from "$module_dir/R/dada2/learnErrors/main.nf" \
    addParams( options: [publish_dir: "denoising"] )
include{ DADA2_DADA } from "$module_dir/R/dada2/dada/main.nf" \
    addParams( options: [publish_dir: "denoising"] )
include{ DADA2_CHIMERA } from "$module_dir/R/dada2/removeBimeraDenovo/main.nf" \
    addParams( options: [publish_dir: "denoising"] )
include{ DADA2_MERGEPAIRS } from "$module_dir/R/dada2/mergePairs/main.nf" \
    addParams( options: [publish_dir: "read-merging"] )

// Other imports
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read-tracking"] )


workflow dada2 {
    take:
    reads // (meta, fwd, rev) tuples

    main:
    tracked_reads = Channel.empty()

    // Remove low quality reads
    qc = DADA2_FILTERANDTRIM( reads )
    derep = DADA2_DEREPFASTQ( qc.fastq ).rds

    // Swith channel to (meta, read)
    derep = derep.transpose().map{ meta, read ->
        def meta_upd = meta.clone() 
        meta_upd['orient'] = (read.getSimpleName() =~ /_R2/ ? "2" : "1")
        [meta_upd, read]
    }
    
    // Remove chimera first for ITS pipeline
    // Requires single end reads or merged paired-end at this stage
    if ( params.early_chimera_removal ) {
        nochim = DADA2_CHIMERA( derep )
        tracked_reads = tracked_reads.mix(nochim.summary)
        derep = nochim.rds
    }

    if ( params.pool != "F" ) {
        derep = derep.map{
            meta, reads ->
            [[id:"all", orient: meta.orient, paired_end: meta.paired_end], reads]
        }.groupTuple(by: 0)
        
        err = DADA2_LEARNERRORS( derep ).rds
        dada = DADA2_DADA( derep.join(err) )
        
    } else {
        // Build Illumina reads error model
        err = DADA2_LEARNERRORS( derep )
        // Denoising
        dada = DADA2_DADA( derep.join(err.rds) )
    }

    // Make raw ASV table
    merged = DADA2_MERGEPAIRS(
        derep.collect{it[1]},
        dada.denoised.collect{it[1]}
    )

    // Read tracking
    merge_summary = SUMMARIZE_TABLE(
        merged.count_table.map{["read-merging", "", it]}
    )
    tracked_reads = tracked_reads.mix(
        qc.summary.mix(dada.summary).mix(merge_summary)
    )

    emit:
    fasta = merged.fasta
    fasta_dup = merged.fasta_dup // only for single_end 
    count_table = merged.count_table
    tracking = tracked_reads
}

