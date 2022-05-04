#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../modules"

// Functions
include { helpMessage; saveParams } from "./util.nf"

// Modules
include{ FASTQC } from "$module_dir/fastqc/main.nf"
include{ GUESS_MATCH_ORDER } from "$module_dir/util/demux/main.nf"
include{ TO_H5 } from "$module_dir/util/demux/main.nf"
include{ ERROR_MODEL } from "$module_dir/util/demux/main.nf"
include{ MAP_INDEX_TO_SAMPLE } from "$module_dir/util/demux/main.nf"
include{ SAMPLE_SIZE_DISTRIBUTION } from "$module_dir/util/demux/main.nf"
include{ WRITE_SAMPLES_TO_FASTQ } from "$module_dir/util/demux/main.nf"
include{ GZIP } from "$module_dir/util/demux/main.nf"


/*
 *
 Beginning of the pipeline
 *
 */

workflow demux {
    take:
    index
    reads
    metadata

    main:

    FASTQC( reads.map{[it[0], it[1..-1]]} )

    barcodes = GUESS_MATCH_ORDER(
        index.map{it[1..-1]},
        metadata
    )

    // Split reads and index for multiprocessing
    index_split = index
        .splitFastq(by: params.n_per_file.toInteger(),
                    file: true,
                    pe: !params.single_barcode)
        .map{[(it[1].getBaseName() =~ /\.([0-9]+)/)[0][1], it[1..-1]]}
    
    reads_split = reads
        .splitFastq(by: params.n_per_file.toInteger(),
                    file: true,
                    pe: !params.single_end)
        .map{[(it[1].getBaseName() =~ /\.([0-9]+)/)[0][1], it[1..-1]]}

    // Convert to h5 for faster loading
    h5 = TO_H5( index_split )

    // Build error model
    model = ERROR_MODEL(
        h5.collect{it[1]},
        barcodes
    )

    // Map index to sample id
    mapping = MAP_INDEX_TO_SAMPLE(
        h5,
        model.h5,
        barcodes
    )

    // Histogram of sample sizes
    SAMPLE_SIZE_DISTRIBUTION( mapping.counts.collect() )

    demux_fq = WRITE_SAMPLES_TO_FASTQ(
        reads_split.join(mapping.tsv)
    )

    GZIP( demux_fq.collect() )
}


workflow {
    index = Channel
        .fromFilePairs("$params.inputdir/*_I{1,2}*.fastq*",
                       size: params.single_barcode ? 1 : 2,
                       flat: true)
        .ifEmpty { error "Cannot find any indexing read in $params.inputdir/" }

    reads = Channel
        .fromFilePairs("${params.inputdir}/*_R{1,2}*.fastq*",
                       size: params.single_end ? 1 : 2,
                       flat: true)
        .ifEmpty { error "Cannot find any sequencing read in $params.inputdir/" }

    metadata = file("$params.inputdir/*.csv")

    demux(index, reads, metadata)
}
