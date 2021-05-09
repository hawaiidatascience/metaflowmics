module_dir = "../modules"

// `preprocess` sub-workflow
include { MOTHUR_ALIGN_SEQS } from "$module_dir/mothur/alignSeqs/main.nf" \
    addParams( options: [publish_dir: "1-MSA-filter"] )
include { MOTHUR_CHIMERA } from "$module_dir/mothur/chimera/main.nf" \
    addParams( options: [publish_dir: "2-chimera-filter"] )
include { MOTHUR_CLASSIFY_SEQS } from "$module_dir/mothur/classifySeqs/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_CLUSTER } from "$module_dir/mothur/cluster/main.nf" \
    addParams( options: [publish_dir: "3-clustering"] )

// `refine` sub-workflow
include { MOTHUR_REMOVE_LINEAGE } from "$module_dir/mothur/removeLineage/main.nf" \
    addParams( options: [publish_dir: "4-lineage-filter"] )
include { MOTHUR_REMOVE_RARE } from "$module_dir/mothur/removeRare/main.nf" \
    addParams( options: [publish_dir: "4-rare-otu-filter"] )
include { GET_SUBSAMPLING_THRESHOLD } from "$module_dir/util/misc/main.nf"
include { MOTHUR_SUBSAMPLE } from "$module_dir/mothur/subsample/main.nf" \
    addParams( options: [publish_dir: "5-subsampling"] )
include { MOTHUR_DIST_SEQS } from "$module_dir/mothur/distSeqs/main.nf" \
    addParams( cutoff: 1-params.lulu_min_match/100, format: "vsearch" )
include { LULU } from "$module_dir/lulu/main.nf" \
    addParams( options: [publish_dir: "6-lulu-filter"] )

// `compile` sub-workflow
include { MOTHUR_CLASSIFY_OTUS } from "$module_dir/mothur/classifyOtus/main.nf"
include { MOTHUR_GET_OTU_REP } from "$module_dir/mothur/getOtuRep/main.nf"
include { MOTHUR_MAKE_DATABASE } from "$module_dir/mothur/makeDatabase/main.nf"

// `sync` sub-workflow
include { MOTHUR_GET_SEQS } from "$module_dir/mothur/getSeqs/main.nf"
include { MOTHUR_GET_OTUS } from "$module_dir/mothur/getOtus/main.nf"

// read tracking
include { MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: "nseqs-sobs" )
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf"


// useful function to split a channel in the different file types
def split_by_extension = branchCriteria {
    fasta: it[-1].getExtension() == "fasta"
    count_table: it[-1].getExtension() == "count_table"
    taxonomy: it[-1].getExtension() == "taxonomy"
    shared: it[-1].getExtension() == "shared"
    summary: it.getExtension() == "summary"
}


workflow preprocess {
    take:
    fasta
    count_table
    db_aln
    db_tax

    main:
    // Filter bad MSA with db
    msa = MOTHUR_ALIGN_SEQS(
        fasta.combine(count_table),
        db_aln
    )

    // Discard chimeric contigs
    chimera = MOTHUR_CHIMERA(
        msa.fasta.combine(msa.count_table)
    )

    // Get contig raw taxonomy
    assignments = MOTHUR_CLASSIFY_SEQS(
        chimera.fasta.combine(chimera.count_table),
        db_aln,
        db_tax
    ).taxonomy

    // Cluster into OTU
    cluster = MOTHUR_CLUSTER(
        chimera.fasta.combine(chimera.count_table),
        params.clustering_thresholds.split(",").collect{it as int}
    )

    // Add OTU identity to channels
    fasta = cluster.list.combine(chimera.fasta).map{[it[0], it[2]]}
    count_table = cluster.list.combine(chimera.count_table).map{[it[0], it[2]]}
    taxonomy = cluster.list.combine(assignments).map{[it[0], it[2]]}

    // File tracking
    tracking = SUMMARIZE_TABLE(msa.count_table.map{["MSA", "", it]}.mix(
        chimera.count_table.map{["chimera-filter", "", it]}
    )).mix(
        MOTHUR_SUMMARY_SINGLE(
            cluster.shared.map{[[step: "clustering", otu_id: it[0]], it[1]]}
        ).summary
    )
    
    emit:
    fasta = fasta
    count_table = count_table
    taxonomy = taxonomy
    list = cluster.list
    shared = cluster.shared
    tracking = tracking
}

workflow refine {
    take:
    count_table
    taxonomy
    list
    repfasta // for LULU
    
    main:
    tracked = Channel.empty()
    
    // Remove unwanted taxa
    if ( params.taxa_to_filter != "" ) {
        (count_table, taxonomy, list, shared) = MOTHUR_REMOVE_LINEAGE(
            count_table.join(taxonomy).join(list)
        )
        tracked = tracked.mix(shared.map{[[step: "taxa-filter", otu_id: it[0]], it[1]]})
    }

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
        dists = MOTHUR_DIST_SEQS( repfasta ).dist
        (shared, repfasta, discard, summary) = LULU(
            dists.join(shared).join(repfasta)
        )
        tracked = tracked.mix(shared.map{[[step: "lulu", otu_id: it[0]], it[1]]})
    }

    tracking = MOTHUR_SUMMARY_SINGLE( tracked ).summary

    emit:
    list = list
    shared = shared
    tracking = tracking
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
