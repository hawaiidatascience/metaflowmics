#!/usr/bin/env nextflow

/*
Backalignment procedure:
 - We collect the alignment at rank r-1
 - ADD_GAPS(): We check with rank r if there are gaps that were deleted --> we add them back
 - COMPUTE_FREQS(): Look at AA frequencies at each position to optimization query/ref mapping
 - UPDATE_MSA: 
*/

nextflow.enable.dsl = 2

process ADD_GAPS {
	label "process_low"
	container "nakor/metaflowmics-python:0.0.1"
	conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

	input:
	path("representative.afa")
	path("representative_aligned.afa")

	output:
	path("representative_aligned_with_gaps.afa")

	script:
	"""
	add_gaps.py --repr representative.afa --upd representative_aligned.afa
	"""
}


process COMPUTE_FREQS {
	label "process_low"
	container "nakor/metaflowmics-python:0.0.1"
	conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

	input:
	path fasta

	output:
	path "freqs*.csv"

	script:
	"""
	#!/usr/bin/env python

	from collections import defaultdict
	from Bio.SeqIO.FastaIO import SimpleFastaParser
	import pandas as pd

	ncols = next(len(x) for (_, x) in SimpleFastaParser(open("$fasta")))
	freqs = [defaultdict(lambda: 0) for _ in range(ncols)]

	with open("$fasta", "r") as reader:
	    for (title, seq) in SimpleFastaParser(reader):
	        for (i, letter) in enumerate(seq):
	            freqs[i][letter] += 1

	freqs = pd.DataFrame(freqs).fillna(0).astype(int)
	freqs.to_csv("freqs_${fasta.getSimpleName()}.csv", index=False)
	"""
}

process UPDATE_MSA {
	tag "$taxa"
	label "process_low"
	publishDir params.outdir, mode: params.publish_dir_mode
	container "nakor/metaflowmics-python:0.0.1"
	conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

	input:
	tuple val(taxa), path(query), path("reference.afa"), path("updated_ref.afa")
	path freqs

	output:
	path "*.updated.afa", emit: afa
	path "inserts.txt", optional: true, emit: txt

	script:
	"""
	update_aln.py \\
	--afa $query \\
	--ref reference.afa \\
	--update updated_ref.afa \\
	--freqs $freqs
	"""
}

workflow update_aln {
	take:
	query
	repr
	repr_aln

	main:

	if (params.field < 3) {
		repr_aln.view{"<${params.field}> $it"}
	}
	
	// Get all reference sequences
	all_repr = repr.collectFile(name: "level-${params.field}.afa").first()
	all_repr_aln = repr_aln.collectFile(name: "level-${params.field}-aln.afa").first()
	
	// Compare references and updated reference to add gaps to the update if they were removed
	all_repr_aln = ADD_GAPS(all_repr, all_repr_aln)
	
	// !!! Then fix UPDATE_MSA to not take that into account
	
	// Get amino acid frequency at each position of the reference
	freqs = COMPUTE_FREQS(all_repr_aln)

	// Get sequences
	repr_per_taxa = all_repr
		.splitFasta(record: [id: true, desc: true, seqString: true])
		.collectFile(){it -> ["${it.desc.tokenize(";")[params.field]}.afa", ">$it.id $it.desc\n$it.seqString\n"]}
		.map{[it.getBaseName(), it]}

	repr_aln_per_taxa = all_repr_aln
		.splitFasta(record: [id: true, desc: true, seqString: true])
		.collectFile(){it -> ["${it.desc.tokenize(";")[params.field]}.afa", ">$it.id $it.desc\n$it.seqString\n"]}
		.map{[it.getBaseName(), it]}

	// Split alignments per taxa
	query_per_taxa = query.map{[it.getName().tokenize(".")[1], it]} // don't forget the field

	// Update the alignment
	query_updated = UPDATE_MSA(
		query_per_taxa.combine(repr_per_taxa, by: 0).combine(repr_aln_per_taxa, by: 0),
		freqs
	)

	query_updated.txt.view()

	emit:
	afa = query_updated.afa
}
