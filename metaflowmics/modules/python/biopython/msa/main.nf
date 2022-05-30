// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process COMPUTE_MSA_REPRESENTATIVE {
    label "process_low"
    publishDir params.outdir, mode: params.publish_dir_mode

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

    input:
    path fasta

    output:
    path "*.repr.fa", emit: repr

    script:
    def prefix = fasta.getSimpleName()
	def skip_gap = params.skip_gap ? "if letter == '-': continue" :""
    """
    #!/usr/bin/env python

    from collections import defaultdict
    import re
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    import pandas as pd

    aln_len = len(next(SimpleFastaParser(open("$fasta")))[1])

    freqs = [defaultdict(lambda: 0) for _ in range(aln_len)]

    with open("$fasta", "r") as reader:
        for (title, seq) in SimpleFastaParser(reader):
            for (i, letter) in enumerate(seq):
	            $skip_gap
	            freqs[i][letter] += 1

    freqs = pd.DataFrame(freqs).fillna(0).astype(int)

    consensus = ''.join(pd.DataFrame(freqs).idxmax(axis=1))
    cons_lineage = title.split()[1]

    with open("${prefix}.repr.fa", "w") as writer:
        writer.write(f">repr|${prefix} {cons_lineage}\\n{consensus}\\n")
    """
}

process UPDATE_MSA_WITH_REF {
    label "process_low"
    publishDir params.outdir, mode: params.publish_dir_mode

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

    input:
    path aln
    path ref

    output:
    path "*.updated.afa", emit: afa

    script:
    """
    #!/usr/bin/env python

    from Bio.SeqIO.FastaIO import SimpleFastaParser

    def map_on_ref(ref, query):
        query_updated = []
        pos = 0
        for aa in ref:
           if aa != '-':
               query_updated.append(query[pos])
               pos += 1
           else:
               query_updated.append('-')
               
        return ''.join(query_updated)

    refs = {}
    with open("$ref") as reader:
        for (title, seq) in SimpleFastaParser(reader):
            lineage = title.split()[1]
            tax = lineage.split("$params.sep")[$params.field]
            refs[tax] = seq

    with open("$aln") as reader:
        for (title, seq) in SimpleFastaParser(reader):
            lineage = title.split()[1]
            tax = lineage.split("$params.sep")[$params.field]
            seq_updated = map_on_ref(refs[tax], seq)
      
            with open(f"{tax}.updated.afa", "a") as writer:
                writer.write(f">{title}\\n{seq_updated}\\n")    
    """    
}

