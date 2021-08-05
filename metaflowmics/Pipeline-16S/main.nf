#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !(params.single_end)
    
module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// DADA2 module imports
include { dada2 } from "$subworkflow_dir/dada2.nf" \
    addParams( outdir: "$params.outdir/interm/read_processing",
              early_chimera_removal: false, format: "mothur" )
// mothur module imports
include { mothur } from "$subworkflow_dir/mothur.nf"

// Other imports
include{ DOWNLOAD_SILVA_FOR_MOTHUR } from "$module_dir/bash/download/main.nf" \
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

    // Read tracking through the pipeline
    READ_TRACKING(
        reads.map{"raw,,${it[0].id},${it[1][0].countFastq()}"}
            .collectFile(newLine: true)
            .mix(asvs.tracking)
            .mix(otus.tracking)
            .collectFile(name: "summary.csv")
    )

    // Visualization
    holoviews(
        otus.shared,
        otus.constaxonomy,
    )

    // Postprocessing
    diversity(
        otus.repfasta,
        otus.shared
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1)
        .map{[ [id: it[0]], it[1] ]}
    saveParams()
    pipeline_16S(reads)
}
