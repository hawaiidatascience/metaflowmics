#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
moduledir = "../../modules"

include{ DOWNLOAD_SILVA_FOR_MOTHUR } from "$moduledir/util/download/main.nf" \
    addParams( db_release: "seed" )
// Mothur module imports
include { MOTHUR_ALIGN_SEQS } from "$moduledir/mothur/alignSeqs/main.nf" \
    addParams( options: [publish_dir: "interm/8-MSA"] )
include { MOTHUR_CHIMERA } from "$moduledir/mothur/chimera/main.nf" \
    addParams( options: [publish_dir: "interm/9-Chimera"] )
include { MOTHUR_CLASSIFY_SEQS } from "$moduledir/mothur/classifySeqs/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_CLUSTER } from "$moduledir/mothur/cluster/main.nf" \
    addParams( options: [publish_dir: "interm/10-Clustering"] )
include { MOTHUR_CLASSIFY_OTUS } from "$moduledir/mothur/classifyOtus/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_GET_SEQS } from "$moduledir/mothur/getSeqs/main.nf"
include { MOTHUR_GET_OTUS } from "$moduledir/mothur/getOtus/main.nf"
include { MOTHUR_MAKE_SHARED } from "$moduledir/mothur/makeShared/main.nf"
include { MOTHUR_SUMMARY_SINGLE } from "$moduledir/mothur/summarySingle/main.nf"
include { MOTHUR_GET_OTU_REP } from "$moduledir/mothur/getOtuRep/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_MAKE_DATABASE } from "$moduledir/mothur/makeDatabase/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_REMOVE_LINEAGE } from "$moduledir/mothur/removeLineage/main.nf" \
    addParams( options: [publish_dir: "interm/11-Filtering"] )
include { MOTHUR_REMOVE_RARE } from "$moduledir/mothur/removeRare/main.nf" \
    addParams( options: [publish_dir: "interm/11-Filtering"] )
include { GET_SUBSAMPLING_THRESHOLD } from "$moduledir/util/misc/main.nf"
include { MOTHUR_SUBSAMPLE } from "$moduledir/mothur/subsample/main.nf" \
    addParams( options: [publish_dir: "interm/12-Subsampling"] )
include { MOTHUR_DIST_SEQS } from "$moduledir/mothur/distSeqs/main.nf" \
    addParams( cutoff: 1-params.lulu_min_match/100, format: 'vsearch' )
// Other imports
include{ SUMMARIZE_TABLE } from "$moduledir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )
include { LULU } from "$moduledir/lulu/main.nf"




workflow mothur_curate {
    take:
    fasta
    count_table

    main:
    // Keep track of the reads in the pipeline
    tracked_ct = Channel.empty()
    
    // Processing with Mothur
    db = DOWNLOAD_SILVA_FOR_MOTHUR()

    // Filter bad MSA with db
    msa = MOTHUR_ALIGN_SEQS(
        fasta.combine(count_table),
        db.align
    )
    tracked_ct = tracked_ct.mix(msa.count_table.map{["MSA", it]})
    
    // Discard chimeric contigs
    chimera = MOTHUR_CHIMERA(
        msa.fasta.combine(msa.count_table)
    )
    tracked_ct = tracked_ct.mix(chimera.count_table.map{["chimera-filter", it]})
    
    // Get contig raw taxonomy
    classify = MOTHUR_CLASSIFY_SEQS(
        chimera.fasta.combine(chimera.count_table),
        db.align,
        db.tax
    )

    fasta = chimera.fasta
    count_table = chimera.count_table
    taxonomy = classify.taxonomy
    
    if ( params.taxa_to_filter != "" ) {
        // Remove unwanted taxa
        lineage = MOTHUR_REMOVE_LINEAGE(
            fasta.concat(count_table, taxonomy).toList()
        )
        fasta = lineage.fasta
        count_table = lineage.count_table
        taxonomy = classify.taxonomy
        tracked_ct = tracked_ct.mix(count_table.map{["taxa-filter", it]})
    }

    // Cluster into OTU
    otus = MOTHUR_CLUSTER(
        fasta.combine(count_table),
        params.clustering_thresholds.split(",").collect{it as int}
    )

    // Discard rare contigs
    rare = MOTHUR_REMOVE_RARE(
        otus.list.combine(count_table)
    )
    
    // Synchronize fasta and taxonomy file
    files = mothur_sync(
        rare.list.combine(fasta.mix(taxonomy)).map{[it[0], it[2]]},
        rare.list
    )

    // Read tracking
    tracked_ct = SUMMARIZE_TABLE( tracked_ct )
    tracked_shared = MOTHUR_SUMMARY_SINGLE(
        otus.shared.map{["clustering-${it[0]}", it[1]]}
            .mix(rare.shared.map{["rare-otus-filter-${it[0]}", it[1]]})
    ).summary
        
    emit:
    fasta=files.fasta
    taxonomy=files.taxonomy
    count_table=rare.count_table
    list=rare.list
    tracking=tracked_ct.mix(tracked_shared)
}

workflow mothur_subsample {
    take:
    count_table
    list

    main:
    subsampling_level = GET_SUBSAMPLING_THRESHOLD(
        count_table.map{it[1]}
    )
    subsampled = MOTHUR_SUBSAMPLE(
        list.join(count_table),
        subsampling_level
    )
    
    tracked_shared = MOTHUR_SUMMARY_SINGLE(
        subsampled.shared.map{["subsampling-${it[0]}", it[1]]}
    )

    emit:
    count_table=subsampled.count_table
    list=subsampled.list
    shared=subsampled.shared
    tracking=tracked_shared.summary
}

workflow mothur_lulu {
    take:
    repfasta
    shared

    main:
    dists = MOTHUR_DIST_SEQS( repfasta ).dist

    lulu = LULU(
        dists.join(shared).join(repfasta)
    )

    emit:
    shared=lulu.abundance
    repfasta=lulu.fasta
    tracking=lulu.summary
}

workflow mothur_sync {
    take:
    files
    list

    main:
    def split_by_extension = branchCriteria {
        count_table: it[1].getExtension() == 'count_table'
        fasta: it[1].getExtension() == 'fasta'
        taxonomy: it[1].getExtension() == 'taxonomy'
    }

    files = MOTHUR_GET_SEQS(
        list.combine(files, by: 0)
    ).branch(split_by_extension)
    
    emit:
    fasta=files.fasta
    count_table=files.count_table
    taxonomy=files.taxonomy
}

workflow mothur_compile {
    take:
    fasta
    count_table
    taxonomy
    list

    main:
    // 1) Abundance table (shared)
    abundance = MOTHUR_MAKE_SHARED(
        list.join(count_table)
    )
    
    // 2) Consensus taxonomy
    cons = MOTHUR_CLASSIFY_OTUS(
        list.join(count_table).join(taxonomy)
    )
    
    // 3) Representative sequences
    rep = MOTHUR_GET_OTU_REP(
        list.join(fasta).join(count_table)
    )
    
    // 4) Database
    if (params.make_db) {
        db = MOTHUR_MAKE_DATABASE(
            list.join(cons.taxonomy).join(rep.fasta).join(rep.count_table)
        ).database
    }
    
    emit:
    shared=abundance.shared
    constaxonomy=cons.taxonomy
    repfasta=rep.fasta
    repcount=rep.count_table    
}

// workflow mothur_postprocessing {
//     take:
    
//     main:
// }
