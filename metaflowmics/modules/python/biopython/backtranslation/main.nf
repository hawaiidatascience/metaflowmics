// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process BACKTRANSLATE {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }
    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::bipython conda-forge::pandas" : null)

    input:
    tuple val(meta), path(afa), path(fna)

    output:
    tuple val(meta), path("*.codons.afa"), emit: fna

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from Bio import SeqIO, AlignIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Data import CodonTable
    from Bio.Align import MultipleSeqAlignment
    from Bio.codonalign import build

    msa = {aln.id: aln for aln in AlignIO.read("$afa", "fasta")}

    dna_all = []
    for dna in SeqIO.parse("$fna", "fasta"):

        if dna.id not in msa:
            continue

        dna.description = msa[dna.id].description

        frame = int(dna.description.split()[-1])
        offset = (frame-1) % 3
        aa_len = (len(dna)-offset) // 3
        dna.seq = dna.seq[offset:3*aa_len+offset]

        if frame > 3:
           dna.seq = dna.seq.reverse_complement()

        dna_all.append(dna)

    msa_clean = MultipleSeqAlignment([
        SeqRecord(
            Seq(str(msa[seq.id].seq.upper()).replace('.', '-')),
            id=seq.id,
            description=seq.description
        ) for seq in dna_all
    ])

    codon_aln = build(msa_clean, dna_all, codon_table=CodonTable.generic_by_id[5])

    with open("${meta.id}.codons.afa", "w") as writer:
        for seq in codon_aln:
            ref = msa[seq.id]
            codons = seq.seq
            seq_str = ''.join(
                '...' if ref.seq[i] == '.'
                else codons.get_codon(i).lower() if ref.seq[i].islower()
                else codons.get_codon(i)
                for i in range(len(ref))
            )
            writer.write(f">{ref.description}\\n{seq_str}\\n")
    """
}
