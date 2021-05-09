#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

// workflows used by mothur main workflow
include { preprocess; refine } from "./mothur-util" \
    addParams( outdir: "$params.outdir/interm/contig_processing" )
include { compile as compile_raw } from "./mothur-util" \
    addParams( outdir: "$params.outdir/raw" )    
include { compile as compile_final } from "./mothur-util" \
     addParams( outdir: "$params.outdir/results" )
include { sync } from "./mothur-util" \
     addParams( outdir: "$params.outdir/results" )


workflow mothur {
    take:
    fasta
    count_table
    db_aln
    db_tax

    main:
    // align.seq, chimera, classify.seq and cluster
    raw = preprocess(
        fasta,
        count_table,
        db_aln,
        db_tax
    )

    // Intermediate results: classify.otus, get.oturep, make.db
    consensus_raw = compile_raw(
        raw.fasta,
        raw.count_table,
        raw.taxonomy,
        raw.list,
        raw.shared
    )

    // remove.lineage, remove.rare, sub.sample and lulu    
    curated = refine(
        raw.count_table,
        raw.taxonomy,
        raw.list,
        consensus_raw.repfasta
    )

    // Subset files to make them coherent
    results = sync(
        raw.fasta.mix(raw.count_table).mix(raw.taxonomy),
        curated.list,
        curated.shared
    )

    // Final results: classify.otus, get.oturep, make.db
    consensus = compile_final(
        results.fasta,
        results.count_table,
        results.taxonomy,
        results.list,
        curated.shared
    )

    // Read tracking
    tracking = raw.tracking.mix(curated.tracking)

    emit:
    repfasta = consensus.repfasta
    constaxonomy = consensus.constaxonomy
    shared = curated.shared
    tracking=tracking
}

