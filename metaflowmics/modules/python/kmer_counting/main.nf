// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process COUNT_KMERS {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process)) }

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::numpy" : null)

    input:
    path fasta

    output:
    path "ref*.npy", emit: kmers
    path "freqs*.npy", emit: freqs

    script:
    if (params.feature == 'nucl') {
        alphabet = "ACGT"
        base = 2
    } else {
        alphabet = "ARNDCQEGHILKMFPSTWYV"
        base = 5
    }
    """
    #!/usr/bin/env python

    from itertools import product
    import numpy as np
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    alphabet = "$alphabet"
    N = len(alphabet)
    k = $params.k

    kmer_idx = {''.join(x): i for i, x in enumerate(product(alphabet, repeat=k))}

    n_db = sum(1 for l in open("$fasta") if l.startswith('>'))

    ref = np.zeros((n_db, N**k), dtype=np.uint16)

    with open("$fasta") as r:
        for i, (_, s) in enumerate(SimpleFastaParser(r)):
            # remove potential gaps
            s = s.replace('-', '')
            # compute kmer counts
            kmer_indices = [kmer_idx[s[i:i+k]] for i in range(len(s)-k+1) if s[i:i+k] in kmer_idx]
            occurrences = np.bincount(kmer_indices, minlength=N**k)

            ref[i] = occurrences

            if i % 1000 == 0:
                print(f"{i:,}/{n_db:,} sequences processed")

    kmer_usage = np.argsort(-(ref > 0).sum(axis=0))

    np.save(f"freqs_{k}mers.npy", kmer_usage)
    np.save(f"ref_{k}mers.npy", ref)
    """
}
