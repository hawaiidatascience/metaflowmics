#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
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
include { TRANSLATE } from "$subworkflow_dir/translation.nf" \
	addParams( outdir: "$params.outdir/translation" )

// Main workflow
workflow download_COI_db {
    take:
    taxa
    
    main:
    // Download iBOL db
    db = DOWNLOAD_IBOL_2(
        taxa
    ).fna.map{it[1]}.collectFile(
		storeDir: params.outdir,
		name: "iBOL_COI_raw.fna"
	)
    
    // Clean missing data in taxonomy and remove duplicated sequences
    db_clean = CLEAN_TAXONOMY(db)

    // Translate sequences
    db_translated = TRANSLATE( db_clean.map{[[id: it.getBaseName()], it]} )

    // Cluster at 98%
    db_no_dup = CDHIT(
        db_translated.faa.map{["db", it]}
    )
}

workflow {
    taxa = Channel.fromPath(params.taxa)
        .splitText()
        .map{it.trim()}
    download_COI_db(taxa)
}
