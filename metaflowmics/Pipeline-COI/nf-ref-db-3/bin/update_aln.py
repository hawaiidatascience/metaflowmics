#!/usr/bin/env python

from pathlib import Path
import argparse
from scipy.special import comb
from functools import lru_cache
from itertools import groupby, combinations

from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--afa', type=str)
    parser.add_argument('--ref', type=str)
    parser.add_argument('--update', type=str)
    parser.add_argument('--freqs', type=str)
    args = parser.parse_args()

    return args

def make_template(consensus, updated, counts):
    ref_itv = get_intervals(consensus)
    other_itv = get_intervals(updated)

    bins = pd.cut(counts.index, pd.IntervalIndex(other_itv))
    counts = counts.groupby(bins).agg(lambda x: list(x)[1:])
    counts.index = pd.IntervalIndex(ref_itv)

    (reference_itv, contracted_idx) = map_intervals(ref_itv, other_itv)

    return (reference_itv, counts, contracted_idx)

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

def map_intervals(ref_itv, other_itv, gap_chars={"-", "."}):
    """
    Map intervals of representative sequences on its aligned version
    Report contracted intervals that need to backpropagate in the other alignments at this level
    """
    # extend contracted intervals
    contracted = other_itv.length < ref_itv.length
    contracted_idx = np.arange(len(other_itv))[contracted]

    # adjust the interval lengths
    ref_itv_adj = pd.Series(ref_itv.length, index=ref_itv)
    ref_itv_adj.loc[~contracted] = other_itv.length[~contracted]

    return (ref_itv_adj, contracted_idx)

def map_on_template(template, query, counts, gap_chars={"-", "."}):
    """
    Map an alignment on the updated template
    The difficulty is that we do not know where new gap were inserted
    among the already existing gaps. There are multiple ways to map.
    We optimize the alignment using the letter frequencies at the 
    undecided positions
    """
    query_aligned = ""

    for (itv, length) in template.iteritems():

        q_sub = query[itv.left:itv.right]
        length_diff = length-len(q_sub)

        if not q_sub:
            query_aligned += '-'*length
            continue

        query_aligned += q_sub[0]
        q_sub = q_sub[1:]            

        # Case 1: template and query have the same length
        if length_diff == 0:
            query_aligned += q_sub
        # Case 2: template is longer than query
        elif length_diff > 0:
            letters = "".join(c for c in q_sub if c not in gap_chars)

            if not letters: # we only have gaps: we just add enough
                query_aligned += "-" * (length-1)
            if letters: # we need to optimize the placement
                opt_aln = find_optimal_mapping(
                    letters, length-1, counts.loc[itv]
                )
                query_aligned += opt_aln
        else: # Should already have been handled when mapping intervals
            raise ValueError("There should not be any contractions at this point")

    return query_aligned

def find_optimal_mapping(letters, total_len, counts):
    print(f"Optimal mapping for n={total_len} and k={len(letters)}")

    counts = np.array([ct for ct in counts.loc[list(letters)]])
    max_score = counts.max(axis=1).sum()
    n_combs = comb(total_len, len(letters), exact=True)

    letter_selector = np.arange(len(letters))
    best_score, best_sol = (0, [])

    for k, positions in enumerate(combinations(range(total_len), len(letters))):
        score = counts[letter_selector, positions].sum()

        if score >= best_score:
            best_score, best_sol = (score, positions)

        if n_combs > 1e6 and (score >= max_score*0.5 or k > 1e7):
            break

    letter_map = dict(zip(positions, letters))
    opt_aln = "".join(letter_map.get(i, '-') for i in range(total_len))

    return opt_aln

if __name__ == '__main__':
    args = parse_args()

    for key in ["ref", "update"]:
        value = getattr(args, key)
        if Path(value).is_file():
            (_, seq) = next(SimpleFastaParser(open(value)))
            setattr(args, key, seq)

    parser = SimpleFastaParser(open(args.afa))

    counts = pd.read_csv(args.freqs)

    (reference_itv, counts, contracted_idx) = make_template(args.ref, args.update, counts)

    output = Path(args.afa).with_suffix(".updated.afa")
    with open(output, "w") as handle:
        for (title, seq) in parser:
            updated_seq = map_on_template(reference_itv, seq, counts)
            handle.write(f">{title}\n{updated_seq}\n")

    if contracted_idx.size > 0:
        with open("inserts.txt", "w") as handle:
            for idx in contracted_idx:
                handle.write(f"{idx}\n")
