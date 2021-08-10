#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../../modules"

// Modules
include{ DOWNLOAD_IBOL } from "$module_dir/bash/download/main.nf" \
    addParams( db_release: params.ibol_release )
include{ CLEAN_TAXONOMY } from "$module_dir/python/taxonomic_cleaner/main.nf"
include{ EMBOSS_TRANSEQ } from "$module_dir/emboss/transeq/main.nf" \
    addParams( table: 5, frames: 6 )
include{ COUNT_KMERS } from "$module_dir/python/kmer_counting/main.nf" \
    addParams( outdir: "$params.outdir/kmer_filtering", k: 3, feature: 'prot', n_sub: 1000 )
include{ KMER_FILTER } from "$module_dir/python/kmer_filter/main.nf" \
    addParams( outdir: "$params.outdir/kmer_filtering", k: 3, feature: 'prot' )
include{ CDHIT } from "$module_dir/cdhit/main.nf" \
    addParams( identity: 1 )

// Main workflow
workflow download_COI_db {
    // Download iBOL db
    versions = Channel.fromList(
        (2..25).collect{String.format('%.2f', it/4)}
    )

    if (!params.test) {
        db = DOWNLOAD_IBOL(versions)
        db_tsv = db.tsv.collectFile(keepHeader: true, skip: 1, storeDir: params.outdir, name: "iBOL_COI.tsv")
    } else {
        db = Channel.fromPath("../../tests/COI/iBOL_COI_sub.tsv")
            .splitCsv(header: true, sep: '\t')
            .map{">$it.seqentryid;k__Animalia,p__$it.phylum_reg,c__$it.class_reg,o__$it.order_reg,f__$it.family_reg,sf__$it.subfamily_reg,g__$it.genus_reg,s__$it.species_reg\n$it.nucraw"}
            .collectFile(newLine: true)
            .map{[0, it]}
            .branch{fna: true}
        db_tsv = Channel.fromPath("../../tests/COI/iBOL_COI_sub.tsv")
    }

    // Clean missing data in taxonomy
    db_clean = CLEAN_TAXONOMY(db.fna)

    // Translate sequences
    db_translated = EMBOSS_TRANSEQ(db_clean)

    // Find frame for genomes with multiple ORF
    db_single = db_translated.single.map{it[1]}.collectFile(name: 'iBOL_COI_single.faa')
    kmer_db = COUNT_KMERS( db_single )
    db_mult = KMER_FILTER(
        db_translated.multiple.map{it[1]}.collectFile(name: 'iBOL_COI_multiple.faa'),
        kmer_db.freqs
    )

    // Combine databases
    db = db_single.mix(db_mult).collectFile(name: 'iBOL_COI_curated.faa', storeDir: params.outdir)

    // Remove duplicates
    db_no_dup = CDHIT(db)
}

workflow {
    download_COI_db()
}
