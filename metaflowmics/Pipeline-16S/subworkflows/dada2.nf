#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !(params.single_end)

moduledir = "../../modules"

include{ DADA2_FILTERANDTRIM } from "$moduledir/dada2/filterAndTrim/main.nf" \
    addParams( options: [publish_dir: "interm/2-Quality_filtering"], trunc_len: 0, trunc_quality: 2, min_read_len: 20 )
include{ DADA2_DEREPFASTQ } from "$moduledir/dada2/derepFastq/main.nf" \
    addParams( options: [publish_dir: "interm/3-Chimera"] )
include{ DADA2_LEARNERRORS } from "$moduledir/dada2/learnErrors/main.nf" \
    addParams( options: [publish_dir: "interm/4-Denoising"] )
include{ DADA2_DADA } from "$moduledir/dada2/dada/main.nf" \
    addParams( options: [publish_dir: "interm/4-Denoising"] )
include{ BUILD_ASV_TABLE } from "$moduledir/util/dada2/main.nf" \
    addParams( options: [publish_dir: "interm/5-Read_Merging"], format: "mothur" )
include{ DADA2_MERGEPAIRS } from "$moduledir/dada2/mergePairs/main.nf" \
    addParams( options: [publish_dir: "interm/5-Read_Merging"] )

// Other imports
include{ SUMMARIZE_TABLE } from "$moduledir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )


workflow dada2 {
    take:
    reads // (meta, fwd, rev) tuples

    main:
    if ( params.single_end ) {
        raw_counts = reads.map{['raw', it[1].countFastq(), it]}
    } else {
        raw_counts = reads.map{['raw', it[1][0].countFastq(), it]}
    }

    // Remove low quality reads
    qc = DADA2_FILTERANDTRIM(
        raw_counts.filter{it[1] > params.min_read_count}.map{it[2]}
    )
    if ( params.single_end ) {
        qc_counts = qc.fastq.map{['QC', it[1].countFastq(), it]}
    } else {
        qc_counts = qc.fastq.map{['QC', it[1][0].countFastq(), it]}
    }
    trimmed_fq = qc_counts.filter{it[1] > params.min_read_count/2}.map{it[2]}
 
    // Build Illumina reads error model
    error_model = DADA2_LEARNERRORS(
        trimmed_fq
    )
    // Denoise reads
    dada = DADA2_DADA(
        trimmed_fq.join(error_model.rds)
    )

    // Make raw ASV table
    if ( params.single_end ) {
        merged = BUILD_ASV_TABLE(
            dada.denoised.collect{it[1]}
        )
    } else {
        merged = DADA2_MERGEPAIRS(
            dada.derep.collect{it[1]},
            dada.denoised.collect{it[1]}
        )
    }

    // Read tracking
    merge_summary = SUMMARIZE_TABLE(
        merged.count_table.map{["read-merging", it]}
    )
            
    tracked_reads = raw_counts.concat(qc_counts)
        .collectFile(name: 'summary.csv', newLine: true){"${it[0]},${it[2][0].id},${it[1]},"}
        .concat(dada.summary, merge_summary)

    emit:
    fasta = merged.fasta
    count_table = merged.count_table
    tracking = tracked_reads
}
