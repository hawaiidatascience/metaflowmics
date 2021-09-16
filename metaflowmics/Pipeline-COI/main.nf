#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]
params.paired_end = !params.single_end

module_dir = "../modules"
subworkflow_dir = "../subworkflows"

// module imports
include{ PEAR } from "$module_dir/pear/main.nf" \
    addParams( options: [publish_dir: "interm/read_processing/merging"] )
include{ CUTADAPT } from "$module_dir/cutadapt/main.nf" \
    addParams( options: [publish_dir: "interm/read_processing/demux"])
include{ READ_TRACKING } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"] )
include{ SYNC_SEQIDS } from "$module_dir/python/biopython/sync_seqids/main.nf"

// mothur imports
include { MOTHUR_ALIGN_SEQS } from "$module_dir/mothur/alignSeqs/main.nf" \
    addParams( options: [publish_dir: "MSA-filter"], no_aln: true )
include { MOTHUR_CHIMERA } from "$module_dir/mothur/chimera/main.nf" \
    addParams( options: [publish_dir: "chimera-filter"] )
include { MOTHUR_CLASSIFY_SEQS } from "$module_dir/mothur/classifySeqs/main.nf" \
    addParams( options: [publish_dir: "raw"] )
include { MOTHUR_CLUSTER } from "$module_dir/mothur/cluster/main.nf" \
    addParams( options: [publish_dir: "clustering"] )
include { MOTHUR_REMOVE_LINEAGE } from "$module_dir/mothur/removeLineage/main.nf" \
    addParams( options: [publish_dir: "lineage-filter"] )
include { MOTHUR_REMOVE_RARE } from "$module_dir/mothur/removeRare/main.nf" \
    addParams( options: [publish_dir: "rare-otu-filter"] )
include { MOTHUR_CLASSIFY_OTUS } from "$module_dir/mothur/classifyOtus/main.nf"
include { MOTHUR_GET_OTU_REP } from "$module_dir/mothur/getOtuRep/main.nf"
include { MOTHUR_SUMMARY_SINGLE } from "$module_dir/mothur/summarySingle/main.nf" \
    addParams( calc: "nseqs-sobs" )
include{ SUMMARIZE_TABLE } from "$module_dir/util/misc/main.nf"

// subworkflow imports
include { dada2 } from "$subworkflow_dir/dada2.nf" \
    addParams( outdir: "$params.outdir/interm/read_processing",
              trunc_len: "0", trunc_quality: 2, min_read_len: 20,              
              early_chimera_removal: false, format: "mothur" )
include { compile as compile_raw } from "$subworkflow_dir/mothur-util" \
    addParams( outdir: "$params.outdir/raw" )
include { compile as compile_final } from "$subworkflow_dir/mothur-util" \
    addParams( outdir: "$params.outdir/results" )
include { translate } from "$subworkflow_dir/translation.nf" \
    addParams( outdir: "$params.outdir/interm/contig_processing/translation",
               db: params.db_aln)
include { hmmer_COI } from "$subworkflow_dir/hmmer.nf" \
    addParams( outdir: "$params.outdir/interm/contig_processing/msa", sync: false )
include { holoviews } from "$subworkflow_dir/holoviews.nf" \
    addParams( options: [publish_dir: "figures"] )
include { diversity } from "$subworkflow_dir/diversity.nf" \
    addParams( options: [publish_dir: "postprocessing"] )

// Main workflow
workflow pipeline_COI {
    take:
    reads
    barcodes
    db_fna
    db_faa
    db_tax

    main:

    // For read tracking through the pipeline
    read_tracking = Channel.empty()
    
    demux = CUTADAPT(
        reads,
        barcodes
    )
    // Convert reads from (dataset, all_reads) to channel emitting (sample_name, reads)
    sample_reads = demux.reads.map{it[1]}.flatten()
        .map{[it.getSimpleName().replaceFirst(/_R[12]$/, ""), it]}
        .groupTuple(by: 0)
        .map{[[id: it[0], paired_end: params.paired_end], it[1]]}

    if (params.single_end || params.merge_with != "PEAR") {
        merged = sample_reads
    }
    else {
        merged = PEAR( sample_reads ).assembled.map{
            {[[id: it[0].id, paired_end: false], it[1]]}
        }
    }

    asvs = dada2( merged )

    // Translation for better alignment
    proteins = translate( asvs.fasta )

    // Remove not translated contigs (i.e. with stop codons)
    synced = SYNC_SEQIDS(
        proteins.faa,
        asvs.fasta.mix(asvs.count_table).collect()
    )

    // alignment with hmmer
    msa = hmmer_COI(
        synced.fna,
        proteins.faa,
        db_faa
    )

    // filter out bad alignments
    msa_filt = MOTHUR_ALIGN_SEQS(
        msa.afa_nucl.combine(synced.count_table),
        db_fna
    )

    // Discard chimeric contigs
    chimera = MOTHUR_CHIMERA(
        msa_filt.fasta.combine(msa_filt.count_table)
    )

    // Taxonomic assignment with Mothur
    mothur_tax = MOTHUR_CLASSIFY_SEQS(
        chimera.fasta.combine(chimera.count_table),
        db_fna,
        db_tax
    ).taxonomy

    // Cluster into OTUs
    cluster = MOTHUR_CLUSTER(
        chimera.fasta.combine(chimera.count_table),
        params.clustering_thresholds.split(",").collect{it as int}
    )

    // Add OTU identity to channels
    list = cluster.list
    shared = cluster.shared
    fasta = list.combine(chimera.fasta).map{[it[0], it[2]]}
    count_table = list.combine(chimera.count_table).map{[it[0], it[2]]}
    taxonomy = list.combine(mothur_tax).map{[it[0], it[2]]}

    interm = compile_raw(
        fasta,
        count_table,
        taxonomy,
        list,
        cluster.shared
    )

    // Remove unwanted taxa
    if ( params.taxa_to_filter != "" ) {
        (count_table, taxonomy, list, shared) = MOTHUR_REMOVE_LINEAGE(
            count_table.join(taxonomy).join(list)
        )
        read_tracking = read_tracking.mix(shared.map{[[step: "taxa-filter", otu_id: it[0]], it[1]]})
    }

    // Discard rare contigs
    (list, shared, count_table) = MOTHUR_REMOVE_RARE(
        list.join(count_table)
    )    
    read_tracking = read_tracking.mix(shared.map{[[step: "rare-otus-filter", otu_id: it[0]], it[1]]})

    // Read tracking through the pipeline
    // READ_TRACKING(
    //     sample_reads.map{"demux,,${it[0].id},${it[1][0].countFastq()}"}
    //         // .mix(merged.map{"pear,,${it[0].id},${it[1].countFastq()}"})
    //         .collectFile(newLine: true)
    //         .mix(asvs.tracking)
    //         .mix(otus.tracking)
    //         .collectFile(name: "summary.csv")
    // )

    // Visualization
    // holoviews(
    //     otus.shared,
    //     otus.constaxonomy,
    // )

    // Postprocessing
    // diversity(
    //     otus.repfasta,
    //     otus.shared
    // )

    
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .map{[ [id: it[0]], it[1] ]}
    barcodes = Channel.fromPath(params.barcodes, checkIfExists: true).collect()

    // Important: DBs need to be value channels 
    db_tax = file(params.db_tax, checkIfExists: true)
    db_fna = file(params.db_aln, checkIfExists: true)
    db_faa = file(params.db_faa, checkIfExists: true)
    
    pipeline_COI(reads, barcodes, db_fna, db_faa, db_tax)
}
