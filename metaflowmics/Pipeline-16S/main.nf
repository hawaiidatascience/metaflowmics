#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.paired_end = !(params.single_end)
params.options = [:]

module_dir = "../modules"
subworkflow_dir = "./subworkflows"

include{ DOWNLOAD_SILVA_FOR_MOTHUR } from "$module_dir/util/download/main.nf" \
    addParams( db_release: params.mothur_db )
include{ READ_TRACKING } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )

// DADA2 module imports
include { dada2 } from "$subworkflow_dir/dada2.nf"
// mothur module imports
include { mothur } from "$subworkflow_dir/mothur.nf"
include { compile } from "$subworkflow_dir/mothur-util.nf"
// Other modules imports
include { holoviews } from "$subworkflow_dir/holoviews.nf"
// Other modules imports
include { diversity } from "$subworkflow_dir/diversity.nf"

// Main workflow
workflow pipeline_16S {
    take:
    reads

    main:
    // Read trimming, QC, denoising and merging
    asvs_ = dada2( reads )
    
    // Download SILVA db for mothur
    db = DOWNLOAD_SILVA_FOR_MOTHUR()

    // OTU curation with mothur
    otus_ = mothur(
        asvs_.fasta,
        asvs_.count_table,
        db.align,
        db.tax
    )

    // Summary
    metagenome_ = compile(
        otus_.fasta,
        otus_.count_table,
        otus_.taxonomy,
        otus_.list,
        otus_.shared
    )

    // Read tracking through the pipeline
    READ_TRACKING(
        asvs_.tracking.mix(otus_.tracking)
            .collectFile(name: 'summary.csv')
    )

    // Visualization
    holoviews(
        otus_.shared,
        metagenome_.constaxonomy,
    )

    // Postprocessing
    diversity(
        metagenome_.repfasta,
        otus_.shared
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1)
        .map{[ [id: it[0], paired: params.paired_end], it[1] ]}
    pipeline_16S(reads)
}
