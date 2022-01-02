// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process CLEAN_TAXONOMY {
    label "process_low"
    publishDir "${params.outdir}/clean_taxonomy", mode: params.publish_dir_mode

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython" : null)

    input:
    path fasta

    output:
    path "*.with_tax.fasta"

    script:
    """
    #!/usr/bin/env python

	import re
	import pandas as pd
	from Bio.SeqIO.FastaIO import SimpleFastaParser

	RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

	def format_taxon_name(taxon):
	    taxon_clean = re.sub("[^A-Za-z0-9_-]", "_", taxon.strip())
	    taxon_clean = re.sub("__+", "_", taxon_clean)
	    return taxon_clean

	def fill_missing_ranks(lineage):
		missing = [i for i,v in enumerate(lineage) if not v][::-1]
		while missing:
			idx = missing.pop()
			prev = lineage[idx-1] # assume the highest level is always known
			prev_rank = RANKS[idx-1][0]
			parts = prev.split("__")

			new_name = f"{prev_rank}__{parts[-1]}"
			if len(parts) > 1:
				prefix = "__".join(parts[:-1])
				new_name = f"{prefix}__{new_name}"
			lineage[idx] = new_name
		return lineage

	criteria = {}
	uniq_seq = {}

	with open("$fasta", "r") as reader:
	    for i, (title, seq) in enumerate(SimpleFastaParser(reader)):
	        info = title.strip().split()
	        lineage = " ".join(info[1:])
	        lineage = lineage.split("$params.sep")
	        lineage = [format_taxon_name(taxon) for taxon in lineage]

	        n_missing = sum(1 for x in lineage if not x)
	        # check if duplicated
	        if seq in uniq_seq and n_missing >= criteria[seq]:
	            continue
	        criteria[seq] = n_missing

	        lineage = "$params.sep".join(
              fill_missing_ranks(lineage)
            )
	        uniq_seq[seq] = f'{info[0]} {lineage}'

	with open("${fasta.getBaseName()}.with_tax.fasta", "w") as writer:
	    for (seq, header) in uniq_seq.items():
	        writer.write(f'>{header}\\n{seq}\\n')
    """
}
