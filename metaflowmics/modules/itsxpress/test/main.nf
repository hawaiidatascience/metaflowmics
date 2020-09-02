#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ITSXPRESS } from '../main.nf' addParams(locus: 'ITS1')

workflow test_single_end {
    def input = []
    input = [
        [ id: 'test', paired: false ],
        [ file("${baseDir}/input/sample1_R1.fastq.gz", checkIfExists: true) ]
    ]

    ITSXPRESS ( input, [ publish_dir:'test_single' ] )
}

workflow test_paired_end {
    def input = []
    input = [
        [ id: 'test', paired: true ],
        [ file("${baseDir}/input/sample1_R1.fastq.gz", checkIfExists: true),
          file("${baseDir}/input/sample1_R2.fastq.gz", checkIfExists: true) ]
    ]

    ITSXPRESS ( input, [ publish_dir:'test_paired' ] )
}

workflow {
    test_single_end()
    test_paired_end()
}
