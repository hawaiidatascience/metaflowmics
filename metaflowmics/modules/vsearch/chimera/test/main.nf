#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_CHIMERA } from '../main.nf'

workflow test {
    def input = []
    input = [
        [ id: 'test' ],
        [ file("${baseDir}/input/dada2_ESVs-100.fasta", checkIfExists: true) ]
    ]

    VSEARCH_CHIMERA ( input, [ publish_dir:'test' ] )
}

workflow {
    test()
}
