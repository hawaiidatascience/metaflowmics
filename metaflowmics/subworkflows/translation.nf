#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.db = null

module_dir = "../modules"


include{ EMBOSS_TRANSEQ } from "$module_dir/emboss/transeq/main.nf" \
    addParams( options: [publish_dir: "transeq"], table: 5, frames: 6 )
include{ COUNT_KMERS } from "$module_dir/python/kmer_counting/main.nf" \
    addParams( options: [publish_dir: "kmer_orf_picking"],
              k: 2, feature: 'prot' )
include{ KMER_FILTER } from "$module_dir/python/kmer_filter/main.nf" \
    addParams( options: [publish_dir: "kmer_orf_picking"],
              k: 2, feature: 'prot', n_sub: 50 )

workflow translate {
    take:
    fasta

    main:
    all_orf = EMBOSS_TRANSEQ( fasta )

    ref = params.db ?
        file(params.db, checkIfExists: true) :
        all_orf.single.map{it[1]}.collectFile(name: "ref.faa")
    
    kmer_db = COUNT_KMERS( ref )

    mult_orf_picked = KMER_FILTER(
        all_orf.multiple,
        kmer_db.freqs
    )

    translated = all_orf.single.mix(mult_orf_picked)

    emit:
    faa = translated
}
