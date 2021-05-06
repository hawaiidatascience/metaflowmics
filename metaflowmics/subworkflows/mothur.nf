#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
module_dir = "../modules"

// workflows used by mothur main workflow
include { sync; compile } from "./mothur-util" \
    addParams( options: [publish_dir: "raw"] )

// mothur module imports
include { MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: "nseqs-sobs" )
include { MOTHUR_ALIGN_SEQS } from "$module_dir/mothur/alignSeqs/main.nf" \
    addParams( options: [publish_dir: "1-MSA-filter"] )
include { MOTHUR_CHIMERA } from "$module_dir/mothur/chimera/main.nf" \
    addParams( options: [publish_dir: "2-chimera-filter"] )
include { MOTHUR_CLASSIFY_SEQS } from "$module_dir/mothur/classifySeqs/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_CLUSTER } from "$module_dir/mothur/cluster/main.nf" \
    addParams( options: [publish_dir: "3-clustering"] )
include { MOTHUR_REMOVE_LINEAGE } from "$module_dir/mothur/removeLineage/main.nf" \
    addParams( options: [publish_dir: "4-lineage-filter"] )
include { MOTHUR_REMOVE_RARE } from "$module_dir/mothur/removeRare/main.nf" \
    addParams( options: [publish_dir: "4-rare-otu-filter"] )

// subsampling
include { GET_SUBSAMPLING_THRESHOLD } from "$module_dir/util/misc/main.nf"
include { MOTHUR_SUBSAMPLE } from "$module_dir/mothur/subsample/main.nf" \
    addParams( options: [publish_dir: "5-subsampling"] )

// lulu
include { MOTHUR_DIST_SEQS } from "$module_dir/mothur/distSeqs/main.nf" \
    addParams( cutoff: 1-params.lulu_min_match/100, format: "vsearch" )
include { LULU } from "$module_dir/lulu/main.nf" \
    addParams( options: [publish_dir: "6-lulu-filter"] )

// Other imports
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf"


// useful function to split a channel in the different file types
def split_by_extension = branchCriteria {
    fasta: it[-1].getExtension() == "fasta"
    count_table: it[-1].getExtension() == "count_table"
    taxonomy: it[-1].getExtension() == "taxonomy"
    shared: it[-1].getExtension() == "shared"
    summary: it.getExtension() == "summary"
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
    (fasta, count_table) = MOTHUR_ALIGN_SEQS(
        fasta.combine(count_table),
        db_aln
    )
    tracked = tracked.mix(count_table.map{["MSA", "", it]})

    // Discard chimeric contigs
    (fasta, count_table) = MOTHUR_CHIMERA(
        fasta.combine(count_table)
    )
    tracked = tracked.mix(count_table.map{["chimera-filter", "", it]})

    // Get contig raw taxonomy
    taxonomy = MOTHUR_CLASSIFY_SEQS(
        fasta.combine(count_table),
        db_aln,
        db_tax
    ).taxonomy

    // Remove unwanted taxa
    if ( params.taxa_to_filter != "" ) {
        (fasta, count_table, taxonomy) = MOTHUR_REMOVE_LINEAGE(
            fasta.concat(count_table, taxonomy).toList()
        )
        tracked = tracked.mix(count_table.map{["taxa-filter", "", it]})
    }

    // Cluster into OTU
    (list, shared) = MOTHUR_CLUSTER(
        fasta.combine(count_table),
        params.clustering_thresholds.split(",").collect{it as int}
    )
    tracked = tracked.mix(shared.map{[[step: "clustering", otu_id: it[0]], it[1]]})

    // Add OTU identity to channels
    fasta = list.combine(fasta).map{[it[0], it[2]]}
    count_table = list.combine(count_table).map{[it[0], it[2]]}
    taxonomy = list.combine(taxonomy).map{[it[0], it[2]]}

    // Intermediate results
    interm = compile(
        fasta,
        count_table,
        taxonomy,
        list,
        shared
    )

    // Discard rare contigs
    (list, shared, count_table) = MOTHUR_REMOVE_RARE(
        list.join(count_table)
    )
    tracked = tracked.mix(shared.map{[[step: "rare-otus-filter", otu_id: it[0]], it[1]]})

    // subsampling
    if (!params.skip_subsampling) {
        subsampling_level = GET_SUBSAMPLING_THRESHOLD( count_table.first().map{it[1]} )
        (shared, list, count_table) = MOTHUR_SUBSAMPLE(
            list.join(count_table),
            subsampling_level
        )
        tracked = tracked.mix(shared.map{[[step: "subsampling", otu_id: it[0]], it[1]]})
    }

    // Co-occurrence pattern correction
    if (!params.skip_lulu) {
        dists = MOTHUR_DIST_SEQS( interm.repfasta ).dist
        (shared, repfasta, discard, summary) = LULU(
            dists.join(shared).join(interm.repfasta)
        )
        tracked = tracked.mix(shared.map{[[step: "lulu", otu_id: it[0]], it[1]]})
    }

    // Subset files to make them coherent
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
    shared=shared
    tracking=tracked_files
}
