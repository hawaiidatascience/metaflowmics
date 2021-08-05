#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DADA2_DADA } from '../main.nf' addParams(truncLen: 200, minLen: 10, truncQ: 2, maxEE: 2, keepPhix: true)

workflow test_fastq {
    def input = []
    input = [
        [ id: 'test' ],
        file("${baseDir}/input/trimmed/sample1_R1.fastq.gz", checkIfExists: true),
        file("${baseDir}/input/errors/sample1_R1-errors.RDS", checkIfExists: true)
    ]

    DADA2_DADA ( input, [ publish_dir:'test_fastq' ] )
}

workflow test_rds {
    def input = []
    input = [
        [ id: 'test' ],
        file("${baseDir}/input/derep/sample1_R1-derep.RDS", checkIfExists: true),
        file("${baseDir}/input/errors/sample1_R1-errors.RDS", checkIfExists: true)
    ]

    DADA2_DADA ( input, [ publish_dir:'test_rds' ] )
}

workflow test_gen {
    reads = Channel.fromFilePairs("${baseDir}/input/trimmed/*_R{1,2}-trimmed.fastq.gz")
    errs = Channel.fromFilePairs("${baseDir}/input/errors/*_R{1,2}-errors.RDS")
    both = reads.combine(errs, by: 0)
        .multiMap {it ->
            fwd: [[id: "${it[0]}_R1"], it[1][0], it[2][0]]
            rev: [[id: "${it[0]}_R2"], it[1][1], it[2][1]]
        }
    input = both.fwd.transpose().mix(both.rev.transpose())
    DADA2_DADA ( input, [ publish_dir:'test' ] )
}

workflow {
    // test_fastq()
    // test_rds()
    test_gen()
}
