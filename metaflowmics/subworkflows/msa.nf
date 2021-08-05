#!/usr/bin/env nextflow

params.options = [:]

moduledir = "../modules"

include{ SPLIT_FASTA } from "$moduledir/python/biopython/split_fasta/main.nf" \
    addParams( sep: params.sep )
include{ COMPUTE_MSA_REPRESENTATIVE ; UPDATE_MSA_WITH_REF } from "$moduledir/python/biopython/msa/main.nf" \
    addParams( sep: params.sep )
include{ MUSCLE } from "$moduledir/muscle/main.nf"



def unknown_filter = branchCriteria {
    unknown: it.getBaseName()[1..-1] == '__'
    known: true
}

workflow msa_agg {
    take:
    fasta

    main:
    // Split fasta in taxonomic groups
    grouped_taxa = SPLIT_FASTA(
        fasta.collectFile(name: "repr_lvl-${params.field}.faa"),
        [field: params.field]
    ).faa.flatten().branch(unknown_filter) // remove unknown taxa
    
    // Align each taxonomic group
    afa_per_taxa = MUSCLE(grouped_taxa.known).afa

    // Get one representative for each aligned fasta to speed up computation
    afa_representatives = COMPUTE_MSA_REPRESENTATIVE(afa_per_taxa).repr

    emit:
    // Re-introduce the sequences for unknown taxa
    afa = afa_per_taxa//.mix(grouped_taxa.unknown)
    repr = afa_representatives.mix(grouped_taxa.unknown)
}

workflow msa_update {
    take:
    level
    fasta
    repr
    
    main:
    // split the repr file to get the lower taxonomic level info
    repr = repr.splitFasta(record: [id: true, seqString: true])
        .map{[it.id.tokenize(';')[-1].tokenize(',')[level], it.seqString]}

    // map it to the corresponding MSA and update the alignment
    fasta = fasta.map{[it.getBaseName(), it]}

    // if (level == 1) {
    //     repr.view{"repr: ${it[0]}"}
    //     fasta.view{"fasta: ${it[0]}"}        
    // }
    
    msa_updated = UPDATE_MSA_WITH_REF(
        fasta.combine(repr, by: 0)
    )

    emit:
    afa = msa_updated
}
