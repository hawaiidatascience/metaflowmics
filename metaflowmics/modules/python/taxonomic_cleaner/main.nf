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

    import re
    import pandas as pd
    from Bio.SeqIO.FastaIO import SimpleFastaParser


    def format_taxon_name(taxon):
        taxon_clean = re.sub("[^A-Za-z0-9_-]", "_", taxon.strip())
        taxon_clean = re.sub("__+", "_", taxon_clean)
        return taxon_clean

    def fill_missing_taxa(lineage):
        cleaned = []
        parent = lineage[0]
        for r, x in enumerate(lineage):
            if not x:
                cleaned.append(f"unknown_{parent}_r{r}")
            else:
                cleaned.append(x)
                parent = x
        return cleaned

    with open("$fasta", "r") as reader, \\
         open("${fasta.getBaseName()}_clean.fasta", "w") as writer:
        for i, (title, seq) in enumerate(SimpleFastaParser(reader)):
            info = title.strip().split()
            lineage = " ".join(info[1:])
            lineage = lineage.split("$params.sep")
            lineage = [format_taxon_name(taxon) for taxon in lineage]
            lineage = "$params.sep".join(
              fill_missing_taxa(lineage)
            )
            writer.write(f'>{info[0]} {lineage}\\n{seq}\\n')
    """
}
