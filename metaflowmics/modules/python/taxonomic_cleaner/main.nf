// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process CLEAN_TAXONOMY {
    label "process_low"
    publishDir "${params.outdir}/clean_taxonomy", mode: params.publish_dir_mode

    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::biopython" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_clean.fasta")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    def format_lineage(lineage):
        cleaned = []
        parent = lineage[0].split('__')[1]
        for x in lineage:
            if x[-2:] == "__":
                cleaned.append(f"{x}unknown_{parent}")
            else:
                cleaned.append(x)
                parent = x.split('__')[1]
        return cleaned

    with open("$fasta", "r") as reader, \\
         open("${fasta.getBaseName()}_clean.fasta", "w") as writer:
        for i, (title, seq) in enumerate(SimpleFastaParser(reader)):
            lineage = ",".join(format_lineage(title.split(",")))
            writer.write(f'>{lineage}\\n{seq}\\n')
    """
}
