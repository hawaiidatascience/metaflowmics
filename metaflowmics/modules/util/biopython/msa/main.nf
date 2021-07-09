// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process COMPUTE_MSA_REPRESENTATIVE {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

    input:
    path fasta

    output:
    path "*.afa", emit: repr

    script:
    """
    #!/usr/bin/env python

    from collections import defaultdict
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    import pandas as pd

    aln_len = len(next(SimpleFastaParser(open("$fasta")))[1])

    freqs = [defaultdict(lambda: 0) for _ in range(aln_len)]

    with open("$fasta", "r") as reader:
        for (title, seq) in SimpleFastaParser(reader):
            for (i, letter) in enumerate(seq):
                if letter != '-':
                    freqs[i][letter] += 1

    consensus = ''.join(pd.DataFrame(freqs).idxmax(axis=1))
    cons_lineage = '$params.sep'.join(title.split(";")[-1].split("$params.sep"))

    with open("${fasta.getBaseName()}.repr.afa", "w") as writer:
        writer.write(f">consensus;{cons_lineage}\\n{consensus}\\n")
    """
}

process UPDATE_MSA_WITH_REF {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge:pandas" : null)

    input:
    tuple val(tax), path(aln), val(ref)

    output:
    path "*.afa", emit: afa

    script:
    prefix = aln.getBaseName()
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

    repr_seq = "$ref"

    with open("$aln") as reader, open("${prefix}_updated.afa", "w") as writer:
        for (title, seq) in SimpleFastaParser(reader):
            seq_updated = map_on_ref(repr_seq, seq)
            writer.write(f">{title}\\n{seq_updated}\\n")
    """    
}
