#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
module_dir = "../../modules"

// mothur module imports
include { MOTHUR_GET_SEQS } from "$module_dir/mothur/getSeqs/main.nf"
include { MOTHUR_GET_OTUS } from "$module_dir/mothur/getOtus/main.nf"
include { MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf"
include { MOTHUR_ALIGN_SEQS } from "$module_dir/mothur/alignSeqs/main.nf" \
    addParams( options: [publish_dir: "interm/8-MSA"] )
include { MOTHUR_CHIMERA } from "$module_dir/mothur/chimera/main.nf" \
    addParams( options: [publish_dir: "interm/9-Chimera"] )
include { MOTHUR_CLASSIFY_SEQS } from "$module_dir/mothur/classifySeqs/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_CLUSTER } from "$module_dir/mothur/cluster/main.nf" \
    addParams( options: [publish_dir: "interm/10-Clustering"] )
include { MOTHUR_CLASSIFY_OTUS } from "$module_dir/mothur/classifyOtus/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_GET_OTU_REP } from "$module_dir/mothur/getOtuRep/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_MAKE_DATABASE } from "$module_dir/mothur/makeDatabase/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_REMOVE_LINEAGE } from "$module_dir/mothur/removeLineage/main.nf" \
    addParams( options: [publish_dir: "interm/11-Filtering"] )
include { MOTHUR_REMOVE_RARE } from "$module_dir/mothur/removeRare/main.nf" \
    addParams( options: [publish_dir: "interm/11-Filtering"] )
include { GET_SUBSAMPLING_THRESHOLD } from "$module_dir/util/misc/main.nf"
include { MOTHUR_SUBSAMPLE } from "$module_dir/mothur/subsample/main.nf" \
    addParams( options: [publish_dir: "interm/12-Subsampling"] )
include { MOTHUR_DIST_SEQS } from "$module_dir/mothur/distSeqs/main.nf" \
    addParams( cutoff: 1-params.lulu_min_match/100, format: 'vsearch' )
// Other imports
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf"
include { LULU } from "$module_dir/lulu/main.nf"



def split_by_extension = branchCriteria {
    fasta: it[-1].getExtension() == 'fasta'
    count_table: it[-1].getExtension() == 'count_table'
    taxonomy: it[-1].getExtension() == 'taxonomy'
    shared: it[-1].getExtension() == 'shared'
    summary: it.getExtension() == 'summary'
}



workflow mothur {
    take:
    fasta
    count_table
    db_aln
    db_tax

    main:
    // Keep track of the reads in the pipeline
    // items are (step, otu_id, file) or a file with those fields
    tracked = Channel.empty()
    
    // Filter bad MSA with db
    msa = MOTHUR_ALIGN_SEQS(
        fasta.combine(count_table),
        db_aln
    )
    tracked = tracked.mix(msa.count_table.map{["MSA", "", it]})
    
    // Discard chimeric contigs
    chimera = MOTHUR_CHIMERA(
        msa.fasta.combine(msa.count_table)
    )
    tracked = tracked.mix(chimera.count_table.map{["chimera-filter", "", it]})
    
    // Get contig raw taxonomy
    classify = MOTHUR_CLASSIFY_SEQS(
        chimera.fasta.combine(chimera.count_table),
        db_aln,
        db_tax
    )

    fasta = chimera.fasta
    count_table = chimera.count_table
    taxonomy = classify.taxonomy
    
    // Remove unwanted taxa
    if ( params.taxa_to_filter != "" ) {
        lineage = MOTHUR_REMOVE_LINEAGE(
            fasta.concat(count_table, taxonomy).toList()
        )
        fasta = lineage.fasta
        count_table = lineage.count_table
        taxonomy = classify.taxonomy
        tracked = tracked.mix(count_table.map{["taxa-filter", "", it]})
    }

    // Cluster into OTU
    otus = MOTHUR_CLUSTER(
        fasta.combine(count_table),
        params.clustering_thresholds.split(",").collect{it as int}
    )
    tracked = tracked.mix(otus.shared.map{["clustering", it[0], it[1]]})

    // Add otu identity to channels
    fasta = otus.list.combine(fasta).map{[it[0], it[2]]}
    count_table = otus.list.combine(count_table).map{[it[0], it[2]]}
    taxonomy = otus.list.combine(taxonomy).map{[it[0], it[2]]}

    // Intermediate results
    interm = compile(
        fasta,
        count_table,
        taxonomy,
        otus.list
    )

    // Discard rare contigs
    rare = MOTHUR_REMOVE_RARE(
        otus.list.join(count_table)
    )
    tracked = tracked.mix(rare.shared.map{["rare-otus-filter", it[0], it[1]]})
    count_table = rare.count_table
    list = rare.list
    shared = rare.shared

    // Subsampling if not skipped
    if (!params.skip_subsampling) {
        subsampled_ = subsample(
            count_table,
            rare.list
        )
        count_table = subsampled_.count_table
        list = subsampled_.list
        shared = subsampled_.shared
        tracked = tracked.mix(subsampled_.tracking)
    }

    // Lulu if not skipped
    if (!params.skip_lulu) {
        lulu_ = lulu(
            interm.repfasta,
            shared
        )
        shared = lulu_.shared
        tracked = tracked.mix(lulu_.tracking)
    }

    files = sync(
        fasta.mix(count_table).mix(taxonomy),
        list,
        shared
    )

    // Read tracking
    tracked_by_ext = tracked.branch(split_by_extension)
    tracked_files = tracked_by_ext.summary.concat(
        SUMMARIZE_TABLE( tracked_by_ext.count_table ),
        MOTHUR_SUMMARY_SINGLE( tracked_by_ext.shared ).summary
    )
   
    emit:
    fasta=files.fasta
    taxonomy=files.taxonomy
    count_table=files.count_table
    list=files.list
    tracking=tracked_files
}

workflow subsample {
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
        subsampled.shared.map{["subsampling", it[0], it[1]]}
    )

    emit:
    count_table=subsampled.count_table
    list=subsampled.list
    shared=subsampled.shared
    tracking=tracked_shared.summary
}

workflow lulu {
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

workflow sync {
    take:
    files
    list
    shared

    main:
    list = MOTHUR_GET_OTUS(
        shared.join(list)
    )

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

    main:
    // 1) Abundance table (shared)
    // abundance = MOTHUR_MAKE_SHARED(
    //     list.join(count_table)
    // )
    
    // 2) Consensus taxonomy
    cons = MOTHUR_CLASSIFY_OTUS(
        list.join(count_table).join(taxonomy)
    )
    
    // 3) Representative sequences
    rep = MOTHUR_GET_OTU_REP(
        list.join(fasta).join(count_table)
    )
    
    // 4) Database
    if (params.compute_mothur_db) {
        db = MOTHUR_MAKE_DATABASE(
            list.join(cons.taxonomy).join(rep.fasta).join(rep.count_table)
        ).database
    }
    
    emit:
    // shared=abundance.shared
    constaxonomy=cons.taxonomy
    repfasta=rep.fasta
    repcount=rep.count_table    
}
