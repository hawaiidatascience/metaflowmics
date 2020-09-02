#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LULU } from '../main.nf' addParams(lulu_min_match: 84, lulu_min_ratio_type: "min", lulu_min_ratio: 1, lulu_min_rel_cooccurence: 1)

workflow test {
    def input = []
    input = [
        100,
        file("${baseDir}/input/matchlist-100.tsv", checkIfExists: true),
        file("${baseDir}/input/dada2_ESVs-100.csv", checkIfExists: true),
    ]

    LULU ( input, [ publish_dir:'test'] )
}


workflow {
    test()
}

