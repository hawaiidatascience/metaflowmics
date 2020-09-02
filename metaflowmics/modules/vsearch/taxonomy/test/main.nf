#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VSEARCH_TAXONOMY } from '../main.nf' addParams(confidence_threshold: 0.5)

workflow test {
    def input = []
    input = [
        [ id:'all_100', otu_id: 100 ],
        file("${baseDir}/input/dada2_ESVs-100.fasta", checkIfExists: true),
        file("${baseDir}/input/dada2_ESVs-100.csv", checkIfExists: true)
    ]
    def db = []
    db = file("${baseDir}/input/test-db/test_unite.fasta")

    VSEARCH_TAXONOMY ( input, db, [ publish_dir:'test' ] )
}

workflow {
    test()
}
