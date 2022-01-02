// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


// Split fasta sequences using the header's taxonomy (split in files corresponding to the same taxonomic rank)

process SPLIT_FASTA {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::biopython" : null)

    input:
    path fasta

    output:
    path "*.main.faa", emit: main
    path "*.others.faa", optional: true, emit: others

    script:
    """
    #!/usr/bin/env python

    import re
    from collections import Counter, defaultdict
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    sequences = defaultdict(lambda: [])
    
    with open("$fasta", "r") as reader:
        for (title, seq) in SimpleFastaParser(reader):
            lineage = title.split()[1]
            tax = lineage.split("$params.sep")[$params.field]
            sequences[tax].append(f">{title}\\n{seq}")

    for (taxa, entries) in sequences.items():        
        seq_str = "\\n".join(entries) + "\\n"

        if len(entries) < $params.min_group_size:
            fname = f"{taxa}.others.faa"
        else:
            fname = f"{taxa}.main.faa"
        
        with open(fname, "w") as writer:
            writer.write(seq_str)
    """
}
