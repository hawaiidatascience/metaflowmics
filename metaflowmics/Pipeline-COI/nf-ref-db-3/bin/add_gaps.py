#!/usr/bin/env python

from pathlib import Path
import argparse
from collections import defaultdict

from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repr', type=str)
    parser.add_argument('--upd', type=str)
    args = parser.parse_args()

    return args

def get_intervals(seq, gap_chars={"-", "."}):
    """
    Get the location and distance between each non-gap characters
    """
    indices = [0]+[i for (i, c) in enumerate(seq) if c not in gap_chars] + [len(seq)]
    intervals = pd.arrays.IntervalArray(
        [pd.Interval(i, j) for (i, j) in zip(indices[:-1], indices[1:])],
        closed='left'
    )
    return intervals

def map_alignment(ref, other, seq_id, gap_chars={"-", "."}):
    """
    Map intervals of representative sequences on its aligned version
    Report contracted intervals that need to backpropagate in the other alignments at this level
    """
    ref_itv = get_intervals(ref)
    other_itv = get_intervals(other)

    # contracted intervals
    contracted = other_itv.length < ref_itv.length
    n_inserts = ref_itv.length - other_itv.length

    inserts = pd.Series(n_inserts[contracted],
                        index=other_itv.right[contracted])

    return inserts

def add_gaps(seq, inserts, seq_id, repres):
    seq_with_gaps = ""
    for i, aa in enumerate(seq):
        if i in inserts:
            seq_with_gaps += "-"*inserts[i]
        seq_with_gaps += aa

    if i+1 in inserts:
        seq_with_gaps += "-"*inserts[i+1]
        
    return seq_with_gaps

if __name__ == '__main__':
    args = parse_args()

    sequences = defaultdict(lambda: {})

    for ftype in ["repr", "upd"]:
        fname = getattr(args, ftype)
        for (title, seq) in SimpleFastaParser(open(fname)):
            sequences[title][ftype] = seq

    inserts = pd.Series(0, index=range(len(seq)+1)) # make sure seq is an updated alignment
    for seq_id, aln in sequences.items():
        new_inserts = map_alignment(aln["repr"], aln["upd"], seq_id).reindex(index=inserts.index)
        more_inserts = inserts < new_inserts
        inserts.loc[more_inserts] = new_inserts[more_inserts]

    inserts = inserts[inserts > 0].astype(int)

    # Update {upd} with new IDs
    with open("representative_aligned_with_gaps.afa", "w") as writer:
        for seq_id, aln in sequences.items():
            with_gaps = add_gaps(aln["upd"], inserts, seq_id, aln["repr"])
            writer.write(f">{seq_id}\n{with_gaps}\n")
