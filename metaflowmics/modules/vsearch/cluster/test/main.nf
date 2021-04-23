#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_CLUSTERING } from '../main.nf'

workflow test {
    def input = []
    input = [
        file("${baseDir}/input/dada2_ESVs-100.fasta", checkIfExists: true)
    ]
    otu_ids = Channel.from(95, 99)
    
    VSEARCH_CLUSTERING ( input, otu_ids, [ publish_dir:'test' ] )
}

workflow {
    test()
}
