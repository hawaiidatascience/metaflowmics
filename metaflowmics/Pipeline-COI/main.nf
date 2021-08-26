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

// subworkflow imports
include { dada2 } from "$subworkflow_dir/dada2.nf" \
    addParams( outdir: "$params.outdir/interm/read_processing",
              trunc_len: "0", trunc_quality: 2, min_read_len: 20,              
              early_chimera_removal: false, format: "mothur" )
include { translate } from "$subworkflow_dir/translation.nf" \
    addParams( outdir: "$params.outdir/interm/contig_processing/translation",
               db: params.db_aln)
include { hmmer_COI } from "$subworkflow_dir/hmmer.nf" \
    addParams( outdir: "$params.outdir/interm/contig_processing/msa" )
include { holoviews } from "$subworkflow_dir/holoviews.nf" \
    addParams( options: [publish_dir: "figures"] )
include { diversity } from "$subworkflow_dir/diversity.nf" \
    addParams( options: [publish_dir: "postprocessing"] )

// Main workflow
workflow pipeline_COI {
    take:
    reads
    barcodes
    db_aln
    // db_tax

    main:
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

    // translation
    asvs_nucl = asvs.fasta.map{[[id: "ASVs"], it]}
    proteins = translate( asvs_nucl )

    // align with hmmer
    alignment = hmmer_COI(
        asvs_nucl,
        proteins.faa,
        db_aln
    )

    // // Translation
    // translated = EMBOSS_TRANSEQ( asvs.fasta )

    // kmer_db = COUNT_KMERS( db_align )
    // translated_mult = KMER_FILTER(
    //     translated.multiple,
    //     kmer_db.freqs
    // )

    // asvs_aa = translated.single.mix(translated_mult)
    //     .collectFile(name: "ASVs.faa")
    //     // .splitFasta(file: true, by: 1000)
    //     .map{[[id: "ASVs"], it]}

    // db_aa = db_align.map{[[id: "db"], it]}
    
    // // Alignment with HMMER
    // hmmdb = HMMER_HMMBUILD( db_align ).hmm

    // hmm_aln_aa = HMMER_HMMALIGN(
    //     asvs_aa.mix(db_aa),
    //     hmmdb
    // )

    // hmm_aln_dna = BACKTRANSLATE(
    //     hmm_aln_aa.afa.join(
    //         asvs.fasta.map{[[id: "ASVs"], it]}.mix(db_dna)
    //     )
    // )
    
    // OTU curation with mothur
    // otus = mothur(
    //     asvs.fasta,
    //     asvs.count_table,
    //     db_align,
    //     db_tax
    // )

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
    align = file(params.db_aln, checkIfExists: true)
    // tax = file(params.db_tax, checkIfExists: true)
    
    pipeline_COI(reads, barcodes, align)
}
