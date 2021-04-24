#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DADA2_MAKESEQUENCETABLE } from '../main.nf'

workflow test {
    def input = []
    input = file("${baseDir}/input/reads_merged-dada2.RDS", checkIfExists: true)

    DADA2_MAKESEQUENCETABLE ( input, [ publish_dir:'test' ] )
}


workflow {
    test()
}
