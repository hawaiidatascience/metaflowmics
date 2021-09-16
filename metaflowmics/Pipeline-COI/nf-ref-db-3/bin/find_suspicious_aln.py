#!/usr/bin/env python

from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import subprocess

import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--aln', type=str)
    parser.add_argument('--outdir', type=str)
    args = parser.parse_args()

    args.aln = Path(args.aln)
    args.outdir = Path(args.outdir)

    return args

args = parse_args()
n_seq = int(subprocess.check_output(["grep", "-c", ">", args.aln]).decode().strip())
n_cols = next(len(x) for (_, x) in SimpleFastaParser(open(args.aln)))

gaps = np.zeros((n_cols, n_seq), dtype=bool)
print("==== Counting gaps ====")
with open(args.aln) as r:
    for j, (title, seq) in enumerate(SimpleFastaParser(r)):
        for i, bp in enumerate(seq):
            gaps[i, j] = bp=="-"
        print(f"{j:,}/{n_seq:,}", end='\r')

sus_all = set()
for (i, col) in enumerate(gaps):
    if (~col).sum() <= 5 or col.sum() <= 5:
        sus = np.arange(n_seq)[~col]
        sus_all.add(sus[0])

print("==== Filtering alignment ====")
outname_filt = Path(args.outdir, f"{args.aln.stem}.kept.faa")
outname_rm = Path(args.outdir, f"{args.aln.stem}.removed.faa")

with open(args.aln) as r, \
     open(outname_filt, "w") as w, \
     open(outname_rm, "w") as rm:
    for j, (title, seq) in enumerate(SimpleFastaParser(r)):
        if j not in sus_all:
            writer = w
        else:
            writer = rm
        writer.write(f">{title}\n{seq.replace('-', '')}\n")
        
        print(f"{j:,}/{n_seq:,}", end='\r')
