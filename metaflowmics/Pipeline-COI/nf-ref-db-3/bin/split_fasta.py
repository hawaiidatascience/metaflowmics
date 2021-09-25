#!/usr/bin/env python

import re
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--field', type=int)    
    parser.add_argument('--min-group-size', type=int, default=2)    
    parser.add_argument('--sep', type=str, default=";")    
    args = parser.parse_args()

    return args

    args = parse_args()

def get_sequences_per_taxa(fasta, field, sep):
    sequences = defaultdict(lambda: [])

    with open(fasta, "r") as reader:
        for (title, seq) in SimpleFastaParser(reader):
            lineage = title.split()[1].split(sep)
            taxa = lineage[field]
            if field > 0:
                parent_taxa = lineage[field-1]
            else:
                parent_taxa = "eukaryota"
            taxa = f"{parent_taxa}.{taxa}"
            sequences[taxa].append(f">{title}\n{seq}")
    return sequences
    
if __name__ == '__main__':
    args = parse_args()

    sequences = get_sequences_per_taxa(args.fasta, args.field, args.sep)

    for (taxa, entries) in sequences.items():
        seq_str = "\n".join(entries) + "\n"

        if len(entries) < args.min_group_size:
            fname = f"{taxa}.others.faa"
        else:
            fname = f"{taxa}.main.faa"

        with open(fname, "w") as writer:
            writer.write(seq_str)
