#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process UPDATE_MSA {
	label "process_low"
	publishDir params.outdir, mode: params.publish_dir_mode
	container "nakor/metaflowmics-python:0.0.1"
	conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

	input:
	tuple path(query), val(ref), val(ref_upd)
	path freqs

	output:
	path "*.updated.afa", emit: afa

	script:
	"""
	update_aln \\
	--afa $query \\
	--repr $ref \\
	--update $ref_upd \\
	--freqs $freqs
	"""
}

process COMPUTE_FREQS {
	label "process_low"
	publishDir params.outdir, mode: params.publish_dir_mode
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
	freqs.to_csv("freqs_${fasta.getBaseName()}.csv")
	"""
}


workflow update_aln {
	take:
	query
	reference
	reference_updated

	main:
	// Get amino acid frequency at each position of the reference
	freqs = COMPUTE_FREQS(reference_updated)
	// Get sequences
	reference_seq = reference
		.splitFasta(record: [id: true, seqString: true ])
		.map{[it.id, it.seqString]}
	reference_upd_seq = reference_updated
		.splitFasta(record: [id: true, seqString: true ])
		.map{[it.id, it.seqString]}
	// Split alignments per taxa
	query_per_taxa = query.map{[it.getName().tokenize(".")[0], it]} // don't forget the field
	// Update the alignment
	reference.view()
	reference_seq.view()
	reference_updated.view()
	reference_upd_seq.view()
	query_per_taxa.view()
	// query.join(reference_seq).join(reference_upd_seq).view()
	query_updated = UPDATE_MSA(
		query_per_taxa.join(reference_seq).join(reference_upd_seq),
		freqs
	)

	emit:
	afa = query_updated
}
