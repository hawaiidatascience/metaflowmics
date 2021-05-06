#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !(params.single_end)
    
module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// DADA2 module imports
include { dada2 } from "$subworkflow_dir/dada2.nf" \
    addParams( options: [publish_dir: "interm/1-read_processing"], early_chimera_removal: false )
// mothur module imports
include { mothur } from "$subworkflow_dir/mothur.nf" \
    addParams( options: [publish_dir: "interm/2-contig_processing"] )
include { compile } from "$subworkflow_dir/mothur-util.nf" \
    addParams( options: [publish_dir: "results"] )
// Other imports
include{ DOWNLOAD_SILVA_FOR_MOTHUR } from "$module_dir/util/download/main.nf" \
    addParams( db_release: params.silva_db )
include{ READ_TRACKING } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )

// Subworkflows
include { holoviews } from "$subworkflow_dir/holoviews.nf" \
    addParams( options: [publish_dir: "figures"] )
include { diversity } from "$subworkflow_dir/diversity.nf" \
    addParams( options: [publish_dir: "postprocessing"] )

// Functions
include { helpMessage ; saveParams } from "./util.nf"


// Main workflow
workflow pipeline_16S {
    take:
    reads

    main:
    // Read trimming, QC, denoising and merging
    asvs = dada2( reads )
    
    // Download SILVA db for mothur
    db = DOWNLOAD_SILVA_FOR_MOTHUR()

    // OTU curation with mothur
    otus = mothur(
        asvs.fasta,
        asvs.count_table,
        db.align,
        db.tax
    )

    // Summary
    metagenome = compile(
        otus.fasta,
        otus.count_table,
        otus.taxonomy,
        otus.list,
        otus.shared
    )

    // Read tracking through the pipeline
    READ_TRACKING(
        reads.map{"raw,,${it[0].id},${it[1].countFastq()}"}
            .collectFile(newLine: true)
            .mix(asvs.tracking)
            .mix(otus.tracking)
            .collectFile(name: "summary.csv")
    )

    // Visualization
    holoviews(
        otus.shared,
        metagenome.constaxonomy,
    )

    // Postprocessing
    diversity(
        metagenome.repfasta,
        otus.shared
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1, flat: true)
        .map{[ [id: it[0], paired: params.paired_end], it[1] ]}
    saveParams()
    pipeline_16S(reads)
}
