#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DADA2_DEREPFASTQ } from '../main.nf'

workflow test {
    def input = []
    input = [
        [ id: 'test' ],
        file("${baseDir}/input/raw/sample1_R1.fastq.gz", checkIfExists: true)
    ]
    input = Channel.fromPath("${baseDir}/input/sample*.fastq.gz")
        .map{[[id: it.getSimpleName()], it]}

    DADA2_DEREPFASTQ ( input, [ publish_dir:'test' ] )
}

workflow {
    test()
}
