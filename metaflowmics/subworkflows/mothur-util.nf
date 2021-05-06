#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
module_dir = "../modules"

// `sync` sub-workflow
include { MOTHUR_GET_SEQS } from "$module_dir/mothur/getSeqs/main.nf"
include { MOTHUR_GET_OTUS } from "$module_dir/mothur/getOtus/main.nf"

// `compile` sub-workflow
include { MOTHUR_CLASSIFY_OTUS } from "$module_dir/mothur/classifyOtus/main.nf"
include { MOTHUR_GET_OTU_REP } from "$module_dir/mothur/getOtuRep/main.nf"
include { MOTHUR_MAKE_DATABASE } from "$module_dir/mothur/makeDatabase/main.nf"

// useful function to split a channel in the different file types
def split_by_extension = branchCriteria {
    fasta: it[-1].getExtension() == "fasta"
    count_table: it[-1].getExtension() == "count_table"
    taxonomy: it[-1].getExtension() == "taxonomy"
    shared: it[-1].getExtension() == "shared"
    summary: it.getExtension() == "summary"
}

workflow sync {
    // Subset all files to make them coherent
    take:
    files
    list
    shared

    main:
    // The shared file is the reference
    // Use it to fix the list file
    list = MOTHUR_GET_OTUS(
        shared.join(list)
    )

    // Fix the other files with the list
    files = MOTHUR_GET_SEQS(
        list.combine(files, by: 0)
    ).branch(split_by_extension)

    emit:
    fasta=files.fasta
    count_table=files.count_table
    taxonomy=files.taxonomy
    list=list
}


workflow compile {
    take:
    fasta
    count_table
    taxonomy
    list
    shared

    main:
    // 1) Consensus taxonomy
    cons = MOTHUR_CLASSIFY_OTUS(
        list.join(count_table).join(taxonomy)
    )

    // 2) Representative sequences
    rep = MOTHUR_GET_OTU_REP(
        list.join(fasta).join(count_table)
    )

    // 3) Database
    if (params.compute_mothur_db) {
        db = MOTHUR_MAKE_DATABASE(
            shared.join(cons.taxonomy).join(rep.fasta).join(rep.count_table)
        ).database
    }

    emit:
    constaxonomy=cons.taxonomy
    repfasta=rep.fasta
    repcount=rep.count_table
}
