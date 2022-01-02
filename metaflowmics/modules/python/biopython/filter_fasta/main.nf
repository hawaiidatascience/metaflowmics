// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process FILTER_FASTA {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython conda-forge::pandas" : null)

    input:
    tuple val(meta), path(fasta), path(hits)

    output:
    path "*.fasta", emit: fasta

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    hits = pd.read_csv("$hits", sep="\\t")
    hits["contig"] = hits.qseqid.str.split("|").str[0]

    kept = set(hits.groupby("contig").pident.agg(lambda x: x.idxmax()))

    with open("$fasta", "r") as reader, \\
         open("${fasta.getBaseName()}_filt.fasta", "w") as writer:
        for i, (title, seq) in enumerate(SimpleFastaParser(reader)):
            if title in kept:
                writer.write(f'{title}\\n{seq}\\n')
    """
}
