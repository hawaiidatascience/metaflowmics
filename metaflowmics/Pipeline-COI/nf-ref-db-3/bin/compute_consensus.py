#!/usr/bin/env python

from collections import defaultdict
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd

import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    args = parser.parse_args()

    return args

def get_consensus(fasta):
    ncols = len(next(SimpleFastaParser(open(args.fasta)))[1])

    freqs = [defaultdict(lambda: 0) for _ in range(ncols)]

    with open(fasta, "r") as reader:
        for (title, seq) in SimpleFastaParser(reader):
            for (i, letter) in enumerate(seq):
                freqs[i][letter] += 1

    freqs = pd.DataFrame(freqs).fillna(0).astype(int)

    consensus = ''.join(pd.DataFrame(freqs).idxmax(axis=1))
    cons_lineage = title.split()[1]

    return (cons_lineage, consensus)


if __name__ == '__main__':
    args = parse_args()

    prefix = Path(args.fasta).with_suffix('').with_suffix('')
    (lineage, seq) = get_consensus(args.fasta)
    
    with open(f"{prefix}.repr.afa", "w") as writer:
        writer.write(f">repr|{prefix} {lineage}\n{seq}\n")

    
