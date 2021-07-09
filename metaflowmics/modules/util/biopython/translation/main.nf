// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process TRANSLATE {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::bipython conda-forge::pandas" : null)

    input:
    path fasta
    path info

    output:
    path "*.faa", emit: faa
    // path "${fasta.getBaseName()}_missing.faa", emit: faa_missing
    // path "summary.csv", emit: summary

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.Seq import Seq
    from Bio.Data import CodonTable

    translation_info = pd.read_csv("$info", sep="\\t", header=None, index_col=0).iloc[:, -1].to_dict()

    def translate(seq, table):
        proteins = []

        for frame in range(0, 3):
            nucl_length = len(seq) - frame
            prot_length = nucl_length // 3
            protein = seq[frame:3*prot_length+frame].translate(table=table)
            if '*' not in protein:
                proteins.append((frame, protein))
        return proteins

    with open("$fasta", "r") as reader, \\
         open("${fasta.getBaseName()}.faa", "w") as writer, \\
         open("${fasta.getBaseName()}_missing.faa", "w") as writer_missing:
        for (title, seq) in SimpleFastaParser(reader):
            lineage = dict(zip(
                ["kingdom", "phylum", "class", "order", "family", "subfamily", "genus", "species"],
                title.split(';')
            ))

            seq = Seq(seq)

            table = translation_info.get(lineage["$params.level"])
            if table is not None:
                proteins = translate(seq, table)
                for (frame, protein) in proteins:
                    writer.write(f">{title} | frame={frame};table={table}\\n{protein}\\n")
            else:
                for table in CodonTable.ambiguous_generic_by_id.keys():
                    proteins = translate(seq, table)
                    for (frame, protein) in proteins:
                        writer_missing.write(f">{title} | frame={frame};table={table}\\n{protein}\\n")
    """
}
