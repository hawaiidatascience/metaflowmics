// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process TRANSLATE {
    tag "$meta"
    label "process_low"
    publishDir params.outdir, mode: params.publish_dir_mode
    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::bipython conda-forge::pandas" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_main.faa"), emit: main
    tuple val(meta), path("*_mult.faa"), emit: mult
    tuple val(meta), path("*_other.faa"), emit: other
    // path "summary.csv", emit: summary

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.Seq import Seq
    from Bio.Data import CodonTable

    def translate(seq, table):
        proteins = []

        for frame in range(0, 3):
            nucl_length = len(seq) - frame
            prot_length = nucl_length // 3
            dna = seq[frame:3*prot_length+frame]
            protein = seq.translate(table=table)
            if '*' not in protein:
                proteins.append((frame, protein))
                continue

            protein = dna.reverse_complement().translate(table=table)
            if '*' not in protein:
                protein.description = "{protein.description} rc"
                proteins.append((frame, protein))
        return proteins

    with open("$fasta", "r") as reader, \\
         open("${fasta.getBaseName()}_main.faa", "w") as writer, \\
         open("${fasta.getBaseName()}_mult.faa", "w") as writer_mult, \\
         open("${fasta.getBaseName()}_other.faa", "w") as writer_other:
        for (title, seq) in SimpleFastaParser(reader):

            if any(x not in "ACGT" for x in seq):
                writer_other.write(f">{title} | unknown_nucleotide \\n{seq}\\n")
                continue

            lineage = dict(zip(
                ["kingdom", "phylum", "class", "order", "family", "subfamily", "genus", "species"],
                title.split(';')
            ))

            seq = Seq(seq)

            proteins = translate(seq, $params.table)

            if len(proteins) == 0:
                proteins = [prot for table in CodonTable.ambiguous_generic_by_id.keys()
                            for prot in translate(seq, table)]
                handle = writer_other
            elif len(proteins) == 1:
                handle = writer
            elif len(proteins) > 1:
                handle = writer_mult
            
            for (frame, protein) in proteins:
                handle.write(f">{title}|frame={frame};table={$params.table}\\n{protein}\\n")
    """
}
