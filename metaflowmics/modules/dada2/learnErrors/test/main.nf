#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DADA2_LEARNERRORS } from '../main.nf'

workflow test_fastq {
    def input = []
    input = [
        [ id: 'test' ],
        file("${baseDir}/input/trimmed/sample1_R1.fastq.gz", checkIfExists: true)
    ]

    DADA2_LEARNERRORS ( input, [ publish_dir:'test_fastq'] )
}

workflow test_rds {
    def input = []
    input = [
        [ id: 'test' ],
        file("${baseDir}/derep/sample1_R1-derep.RDS", checkIfExists: true)
    ]

    DADA2_LEARNERRORS ( input, [ publish_dir:'test_rds'] )
}


workflow test_gen {
    input = Channel.fromFilePairs("${baseDir}/input/trimmed/sample*_R{1,2}*.fastq.gz")
        .multiMap{it ->
            fwd: [[id: "${it[0]}_R1"], it[1][0]]
            rev: [[id: "${it[0]}_R2"], it[1][1]]
        }
    DADA2_LEARNERRORS ( input.fwd.mix(input.rev), [ publish_dir:'test' ] )
}


workflow {
    // test_rds()
    // test_fastq()
    test_gen()
}

