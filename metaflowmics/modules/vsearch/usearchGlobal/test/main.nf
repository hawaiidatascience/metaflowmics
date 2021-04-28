#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_USEARCHGLOBAL } from '../main.nf' addParams(lulu_min_match: 0.84)

workflow test {
    def input = [100, file("${baseDir}/input/dada2_ESVs-100.fasta", checkIfExists: true)]
    
    VSEARCH_USEARCHGLOBAL ( input, [ publish_dir:'test' ] )
}

workflow {
    test()
}
