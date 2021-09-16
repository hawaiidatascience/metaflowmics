#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.options = [:]
moduledir = "../modules"

include{ SPLIT_FASTA } from "$moduledir/python/biopython/split_fasta/main.nf"
include{ MUSCLE } from "$moduledir/muscle/main.nf"
include{ COMPUTE_MSA_REPRESENTATIVE } from "$moduledir/python/biopython/msa/main.nf" \
    addParams( outdir: "$params.outdir/taxa_representative" )


// Main workflow
workflow align_tax_groups {
    take:
    fasta

    main:
    // Split fasta in taxonomic groups
    grouped_taxa = SPLIT_FASTA(
        fasta.collectFile(name: "repr_lvl-${params.field}.faa")
    )

    // Align each taxonomic group
    afa_per_taxa = MUSCLE(grouped_taxa.main.flatten()).afa

    // Get one representative for each aligned fasta to speed up computation
    afa_representatives = COMPUTE_MSA_REPRESENTATIVE(afa_per_taxa).repr

    emit:
    // Re-introduce the sequences for taxa with one sequence
    afa = afa_per_taxa.mix(grouped_taxa.others.flatten())
    repr = afa_representatives.mix(grouped_taxa.others.flatten())
}

