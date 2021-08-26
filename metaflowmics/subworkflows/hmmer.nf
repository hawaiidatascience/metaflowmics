#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../modules"

include{ HMMER_HMMBUILD } from "$module_dir/hmmer/hmmbuild/main.nf" \
    addParams( options: [args: "--amino"] )
include{ HMMER_HMMALIGN } from "$module_dir/hmmer/hmmalign/main.nf" \
    addParams( options: [publish_dir: "hmmalign", args: "--amino"] )
include{ BACKTRANSLATE } from "$module_dir/python/biopython/backtranslation/main.nf" \
    addParams( options: [publish_dir: "backtranslate"] )


workflow hmmer_COI {
    take:
    fna // unaligned DNA sequence
    faa // corresponding protein sequence
    db // protein reference MSA

    main:
    // Build HMM model
    hmmdb = HMMER_HMMBUILD( db ).hmm

    // Align proteins against database
    msa_prot = HMMER_HMMALIGN(
        faa,
        hmmdb
    )

    // Reverse translate to get the DNA alignment
    msa_nucl = BACKTRANSLATE(
        msa_prot.afa.join(fna)
    )
    
    emit:
    afa_nucl = msa_nucl.fna
    afa_prot = msa_prot.afa
}

workflow {
    hmmer_COI(
        Channel.fromPath(params.fna).map{[[id: "db"], it]},
        Channel.fromPath(params.faa).map{[[id: "db"], it]},
        Channel.fromPath(params.db)
    )
}
