#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.paired_end = !(params.single_end)
params.options = [:]

moduledir = "../../modules"

include { compile } from "./mothur.nf"
include{ READ_TRACKING } from "$moduledir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )


workflow summarize {
    take:
    fasta
    count_table
    taxonomy
    list
    tracking

    main:

    // consensus files from ASV to OTU level
    metagenome = compile(
        fasta,
        count_table,
        taxonomy,
        list
    )

    // Read tracking through the pipleine
    summary = READ_TRACKING( tracking.collectFile(name: 'summary.csv') )    
    
    // Summary plots
    
}
