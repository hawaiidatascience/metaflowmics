#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.paired_end = !(params.single_end)
params.options = [:]

// DADA2 module imports
include { dada2 } from './subworkflows/dada2.nf'
// Mothur module imports
include { mothur_curate; mothur_compile; mothur_subsample; mothur_lulu; mothur_sync } \
    from './subworkflows/mothur.nf'
// Other imports
include{ MOTHUR_GET_OTUS } from '../modules/mothur/getOtus/main.nf' \
    addParams( options: [publish_dir: 'read_tracking'] )
include{ READ_TRACKING } from '../modules/util/misc/main.nf' \
    addParams( options: [publish_dir: 'read_tracking'] )

// Workflows definition
workflow pipeline_16S {
    take:
    reads

    main:
    // Read trimming, QC, denoising and merging
    asvs_ = dada2( reads )
    
    // Filtering of OTU (multiple criteria)
    otus_ = mothur_curate(
        asvs_.fasta,
        asvs_.count_table
    )

    // Intermediate results
    interm = mothur_compile(
        otus_.fasta,
        otus_.count_table,
        otus_.taxonomy,
        otus_.list
    )

    // Keep track of the following files as we go through optional steps
    count_table = otus_.count_table
    shared = interm.shared
    tracking = asvs_.tracking.mix(otus_.tracking)

    // Subsampling if not skipped
    if (!params.skip_subsampling) {
        subsampled_ = mothur_subsample(
            count_table,
            otus_.list
        )
        count_table = subsampled_.count_table
        shared = subsampled_.shared
        tracking = tracking.mix(subsampled_.tracking)
    }

    // Lulu if not skipped
    if (!params.skip_lulu) {
        lulu_ = mothur_lulu(
            interm.repfasta,
            shared
        )
        shared = lulu_.shared
        tracking = tracking.mix(lulu_.tracking)
    }

    // Sync all files
    list = MOTHUR_GET_OTUS(
        shared.join(otus_.list)
    )
    files = mothur_sync(
        otus_.fasta.mix(count_table).mix(otus_.taxonomy),
        list
    )

    summary = READ_TRACKING( tracking.collectFile() )    
}


workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1)
        .map{[ [id: it[0], paired: params.paired_end], it[1] ]}
    pipeline_16S(reads)
}
