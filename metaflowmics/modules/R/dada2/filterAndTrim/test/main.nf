#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DADA2_FILTERANDTRIM } from '../main.nf' addParams(truncLen: "250, 170", minLen: 20, truncQ: 2, maxEE: 3, keepPhix: true)

workflow test_single_end {
    def input = []
    input = [
        [ id: 'test', paired: false ],
        [ file("${baseDir}/input/raw/sample1_R1.fastq.gz", checkIfExists: true) ]
    ]

    DADA2_FILTERANDTRIM ( input, [ publish_dir:'test_single' ] )
}

workflow test_paired_end {
    def input = []
    input = [
        [ id: 'test', paired: true ],
        [ file("${baseDir}/input/raw/sample1_R1.fastq.gz", checkIfExists: true),
         file("${baseDir}/input/raw/sample1_R2.fastq.gz", checkIfExists: true)]
    ]

    DADA2_FILTERANDTRIM ( input, [ publish_dir:'test_paired' ] )
}

workflow test_gen {
    input = Channel.fromFilePairs("${baseDir}/input/raw/*_R{1,2}.fastq.gz")
        .map{[[id: it[0], paired: true], it[1]]}
    DADA2_FILTERANDTRIM ( input, [ publish_dir:'test_gen' ] )
}

workflow {
    // test_single_end()
    // test_paired_end()
    test_gen()
}
