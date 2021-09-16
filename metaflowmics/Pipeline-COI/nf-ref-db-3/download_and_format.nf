#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.custom_db = ""

module_dir = "../../modules"
subworkflow_dir = "../../subworkflows"

// Modules
include{ DOWNLOAD_IBOL_2 } from "$module_dir/bash/download/main.nf" \
    addParams( db_release: params.ibol_release )
include{ CLEAN_TAXONOMY } from "$module_dir/python/taxonomic_cleaner/main.nf"
include{ CDHIT } from "$module_dir/cdhit/main.nf" \
    addParams( identity: 0.98 )

include { translate } from "$subworkflow_dir/translation.nf"


// Main workflow
workflow download_COI_db {
    take:
    taxa
    
    main:
    // Download iBOL db
    db = DOWNLOAD_IBOL_2(
        taxa
    )
    
    // Clean missing data in taxonomy
    db_clean = CLEAN_TAXONOMY(db.fna).map{it[1]}

    // Save the fna to a single file for later use
    db_clean.collectFile(
        storeDir: params.outdir,
        name: "iBOL_COI.fna"
    )

    // Translate sequences
    db_translated = translate( db_clean )

    // Remove duplicates
    db_no_dup = CDHIT(
        db_translated.faa
            .collectFile(name: "iBOL_COI.faa")
    )
}

workflow {
    taxa = Channel.fromPath(params.taxa)
        .splitText()
        .map{it.trim()}
    download_COI_db(taxa)
}
