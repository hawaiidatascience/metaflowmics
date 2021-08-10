// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process KMER_FILTER {
    label "process_high"
    publishDir params.outdir, mode: "copy"
    
    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::numpy" : null)

    input:
    path query
    path db

    output:
    path "*.fasta"

    script:
    prefix = "${query.getBaseName()}_curated"

    if (params.feature == 'nucl') {
        alphabet = "ACGT"
        base = 2
    } else {
        alphabet = "ARNDCQEGHILKMFPSTWYV"
        base = 5
    }
    """
    #!/usr/bin/env python

    from itertools import product, groupby
    import numpy as np
    from Bio.SeqIO.FastaIO import SimpleFastaParser


    alphabet = "$alphabet"
    N = len(alphabet)
    k = $params.k

    kmer_usage = set(np.load("$db")[:100])
    kmer_db = {''.join(x) for i, x in enumerate(product(alphabet, repeat=k)) if i in kmer_usage}

    with open("$query") as r, open("${prefix}.fasta", "w") as w:
        for i, (_, entries) in enumerate(
            groupby(SimpleFastaParser(r), lambda x: x[0].split()[0])
        ):
            best_score = 0
            selected = ()

            for (t, s) in entries:
                # compute kmer counts
                score = sum(1 for i in range(len(s)-k+1) if s[i:i+k] in kmer_db)

                if score > best_score:
                    selected = (t.split('|')[0], s)
                    best_score = score

            w.write('>{}\\n{}\\n'.format(*selected))

            if i % 1000 == 0:
                print(f"{i:,} sequences processed")
    """
}