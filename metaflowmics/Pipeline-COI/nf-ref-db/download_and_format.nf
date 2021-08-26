#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../../modules"
subworkflow_dir = "../../subworkflows"

// Modules
include{ DOWNLOAD_IBOL } from "$module_dir/bash/download/main.nf" \
    addParams( db_release: params.ibol_release )
include{ CLEAN_TAXONOMY } from "$module_dir/python/taxonomic_cleaner/main.nf"
include{ CDHIT } from "$module_dir/cdhit/main.nf" \
    addParams( identity: 0.98 )

include { translate } from "$subworkflow_dir/translation.nf"


// Main workflow
workflow download_COI_db {    
    // Download iBOL db
    versions = Channel.fromList(
        (2..25).collect{String.format('%.2f', it/4)}
    )

    db = DOWNLOAD_IBOL(versions)
    db_tsv = db.tsv.collectFile(
        keepHeader: true,
        skip: 1,
        storeDir: "$params.outdir/download",
        name: "iBOL_COI.tsv"
    )
    
    // Clean missing data in taxonomy
    db_clean = CLEAN_TAXONOMY(db.fna)

    // Translate sequences
    db_translated = translate( db_clean )

    // Remove duplicates
    db_no_dup = CDHIT(
        db_translated.faa
            .map{it[1]}
            .collectFile(name: "iBOL_COI.faa")
            .map{[[id: "db"], it]}
    )
}

workflow {
    download_COI_db()
}
