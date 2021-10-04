#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.sync = true

module_dir = "../modules"

include{ HMMER_HMMBUILD } from "$module_dir/hmmer/hmmbuild/main.nf" \
    addParams( options: [args: "--amino"] )
include{ HMMER_HMMALIGN } from "$module_dir/hmmer/hmmalign/main.nf" \
    addParams( options: [args: "--amino", publish_dir: "hmmalign"] )
include{ EMBOSS_TRANALIGN } from "$module_dir/emboss/tranalign/main.nf" \
    addParams( options: [publish_dir: "tranalign"] )
include{ SYNC_SEQIDS } from "$module_dir/python/biopython/sync_seqids/main.nf"


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
    ).afa

    if (params.sync) {
        fna = SYNC_SEQIDS(
            msa_prot,
            fna.collect{it[1]} // make it a list for SYNC_SEQIDS
        ).fna
    }
    
    // Reverse translate to get the DNA alignment
    msa_nucl = EMBOSS_TRANALIGN(
        msa_prot.join(fna)
    )
    
    emit:
    afa_nucl = msa_nucl.fna
    afa_prot = msa_prot
}

workflow {
    hmmer_COI(
        Channel.fromPath(params.fna),
        Channel.fromPath(params.faa),
        file(params.db, checkIfExists: true)
    )
}
