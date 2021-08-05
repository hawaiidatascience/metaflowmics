#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// module imports
include{ PEAR } from "$module_dir/pear/main.nf" \
    addParams( options: [publish_dir: "interm/read_merging"] )
include{ CUTADAPT } from "$module_dir/cutadapt/main.nf" \
    addParams( options: [publish_dir: "interm/cutadapt"], single_end: true )

// subworkflow imports

// Main workflow
workflow pipeline_COI {
    take:
    reads
    barcodes

    main:
    merged = PEAR( reads )
    demux = CUTADAPT(
        merged.assembled,
        barcodes
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .map{[ [id: it[0]], it[1] ]}
    barcodes = file(params.barcodes, checkIfExists: true)

    pipeline_COI(reads, barcodes)
}
