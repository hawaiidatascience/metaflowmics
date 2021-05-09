#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.options = [:]

module_dir = "../modules"
subworkflow_dir = "../subworkflows"

include{ DOWNLOAD_UNITE } from "$module_dir/util/download/main.nf" \
    addParams( db_release: "fungi" )
include{ ITSXPRESS } from "$module_dir/itsxpress/main.nf" \
    addParams( options: [publish_dir: "interm/read_processing/itsxpress"] )
include{ VSEARCH_CLUSTER } from "$module_dir/vsearch/cluster/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/clustering"] )
include{ VSEARCH_USEARCH_GLOBAL } from "$module_dir/vsearch/usearchGlobal/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/lulu"] )
include{ LULU } from "$module_dir/lulu/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/lulu"] )
include{ DADA2_ASSIGN_TAXONOMY } from "$module_dir/dada2/assignTaxonomy/main.nf" \
    addParams( options: [publish_dir: "interm/contig_processing/taxonomy"] )
include{ SUMMARIZE_TABLE; READ_TRACKING } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "read_tracking"], taxa_are_rows: "T" )
include{ CONVERT_TO_MOTHUR_FORMAT } from "$module_dir/util/misc/main.nf" \
    addParams( options: [publish_dir: "results"], taxa_are_rows: "T" )

// subworkflow imports
include { dada2 } from "$subworkflow_dir/dada2.nf" \
    addParams( outdir: "$params.outdir/interm/read_processing",
              trunc_len: 0, trunc_quality: 2, min_read_len: 20, paired_end: false,
              early_chimera_removal: true, format: "VSEARCH" )
include { holoviews } from "$subworkflow_dir/holoviews.nf" \
    addParams( options: [publish_dir: "figures"] )
include { diversity } from "$subworkflow_dir/diversity.nf" \
    addParams( options: [publish_dir: "postprocessing"], skip_unifrac: true, unifrac: '' )

// Functions
include { helpMessage; saveParams } from "./util.nf"


// Main workflow
workflow pipeline_ITS {
    take:
    reads

    main:
    tracked_files = Channel.empty()
    raw_counts = reads.map{["raw", it[1][0].countFastq(), it]}

    // ===== Extract ITS marker =====
    its = ITSXPRESS(
        raw_counts.filter{it[1] > params.min_read_count}.map{it[2]}
    )

    // ===== Denoising into ASVs =====
    asvs = dada2(its.fastq)
    
    // ===== Cluster ASVs into OTUs =====
    otus = VSEARCH_CLUSTER(
        asvs.fasta_dup,
        params.clustering_thresholds.split(",").findAll({it!="100"}).collect{it as int}
    )
    otus_fasta = asvs.fasta.map{[100, it]}.mix(otus.fasta)
    otus_table = asvs.count_table.map{[100, it]}.mix(otus.tsv)

    otus_summary = SUMMARIZE_TABLE(otus_table.map{["clustering", it[0], it[1]]})

    // ===== Optional co-occurrence pattern correction =====
    if (!params.skip_lulu) {
        matchlist = VSEARCH_USEARCH_GLOBAL(
            otus_fasta
        )
        lulu = LULU(
            matchlist.tsv.join(otus_table).join(otus_fasta)
        )
        tracked_files = tracked_files.mix(lulu.summary)
        otu_fastas = lulu.fasta
        otu_tables = lulu.abundance
    }

    // ===== Taxonomy assignment with Sintax =====
    unite_db = DOWNLOAD_UNITE()
    taxonomy = DADA2_ASSIGN_TAXONOMY(
        otus_fasta,
        unite_db
    ).taxonomy

    // ===== Track read in the pipeline =====    
    tracked_files = tracked_files
        .mix(
            raw_counts.map{"raw,,${it[2][0].id},${it[1]}"}
                .mix(its.fastq.map{"itsxpress,,${it[0].id},${it[1].countFastq()}"})
                .collectFile(newLine: true)
        )
        .mix(asvs.tracking)
        .mix(otus_summary)
        .collectFile(name: "summary.csv")

    summary = READ_TRACKING( tracked_files )

    // Conversion to mothur format
    mothur_files = CONVERT_TO_MOTHUR_FORMAT(
        otus_table.join(taxonomy)
    )

    // ===== Diversity =====
    diversity(
        otus_fasta,
        mothur_files.shared
    )
    
    // ===== Plots =====
    holoviews(
        mothur_files.shared,
        mothur_files.taxonomy
    )
}

workflow {
    reads = Channel.fromFilePairs(params.reads, size: params.paired_end ? 2 : 1)
        .map{[ [id: it[0]], it[1] ]}

    saveParams()
    pipeline_ITS(reads)
}
