#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DADA2_MERGEPAIRS } from '../main.nf' addParams(min_overlap: 10, max_mismatch: 1)

workflow test {
    derep = Channel.fromPath("${baseDir}/input/derep/sample*-derep.RDS")
    denoised = Channel.fromPath("${baseDir}/input/denoised/sample*-denoised.RDS")

    DADA2_MERGEPAIRS ( derep.collect(), denoised.collect(), [ publish_dir:'test' ] )
}

workflow {
    test()
}
