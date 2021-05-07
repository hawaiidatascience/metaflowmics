#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.early_chimera_removal = params.early_chimera_removal ?: false

module_dir = "../modules"

include{ DADA2_FILTERANDTRIM } from "$module_dir/dada2/filterAndTrim/main.nf" \
    addParams( options: [publish_dir: "1-quality-filtering"] )
include{ DADA2_LEARNERRORS } from "$module_dir/dada2/learnErrors/main.nf" \
    addParams( options: [publish_dir: "2-denoising"] )
include{ DADA2_DADA } from "$module_dir/dada2/dada/main.nf" \
    addParams( options: [publish_dir: "2-denoising"] )
include{ BUILD_ASV_TABLE } from "$module_dir/util/dada2/main.nf" \
    addParams( options: [publish_dir: "3-read-merging"] )
include{ DADA2_MERGEPAIRS } from "$module_dir/dada2/mergePairs/main.nf" \
    addParams( options: [publish_dir: "3-read-merging"] )

// If VSEARCH is used for chinera
include{ DADA2_DEREPFASTQ } from "$module_dir/dada2/derepFastq/main.nf" \
    addParams( options: [publish_dir: "2.5-Chimera"] )
include{ VSEARCH_CHIMERA } from "$module_dir/vsearch/chimera/main.nf" \
    addParams( options: [publish_dir: "2.5-Chimera"] )
include{ SUBSET_READS_RDS } from "$module_dir/util/dada2/main.nf" \
    addParams( options: [publish_dir: "2.5-OTU_clustering"] )

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

    // Dereplicate reads
    derep = DADA2_DEREPFASTQ( qc.fastq )
    
    // Remove chimera first for ITS pipeline
    // Ask Anthony if dada2 removeBimeraDenovo is OK
    if ( params.early_chimera_removal ) {
        chimera_filt = VSEARCH_CHIMERA( derep.fasta )
        chimera_filt = SUBSET_READS_RDS( derep.rds.join(chimera_filt.fasta) )
        tracked_reads = tracked_reads.mix(chimera_filt.summary)
    } else {
        chimera_filt = derep
    }

    // Build Illumina reads error model
    error_model = DADA2_LEARNERRORS( chimera_filt.rds )
    // Denoising
    dada = DADA2_DADA( chimera_filt.rds.join(error_model.rds) )

    // Make raw ASV table
    if ( params.paired_end ) {
        merged = DADA2_MERGEPAIRS(
            chimera_filt.rds.collect{it[1]},
            dada.denoised.collect{it[1]}
        )
        fasta_dup = Channel.empty()
    } else {
        merged = BUILD_ASV_TABLE(
            dada.denoised.collect{it[1]}
        )
        fasta_dup = merged.fasta_dup
    }

    // Read tracking
    merge_summary = SUMMARIZE_TABLE(
        merged.count_table.map{["read-merging", "", it]}
    )
    tracked_reads = tracked_reads.mix(
        qc.summary.mix(dada.summary).mix(merge_summary)
    )

    emit:
    fasta = merged.fasta
    fasta_dup = fasta_dup // only for single_end 
    count_table = merged.count_table
    tracking = tracked_reads
}
