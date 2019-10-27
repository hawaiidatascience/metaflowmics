#!/usr/bin/env python3

'''
- Read index and barcodes
- Match each read_id to a sample_id
- Write the result in a file
'''

import argparse
from itertools import product

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

NUCLS = list('ACGTN')
MAPPING = {ord(nucl): i for i, nucl in enumerate(NUCLS)}
MAPPING_BASE5 = str.maketrans(''.join(NUCLS), '01234')

def parse_args():
    '''
    Get command line arguments
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--fwd', type=str, help='Path to fwd index')
    parser.add_argument('--rev', type=str, help='Path to rev index')
    parser.add_argument('--meta', type=str, help='Path to metadata')
    parser.add_argument('--max_mismatches', type=int, default=1)
    parser.add_argument('--revMapping', action='store_true', default='false')
    parser.add_argument('--strategy', type=str, choices=['best', 'min_all', 'discard'], default='best', help='Strategy to handle multimappers')    
    args = parser.parse_args()

    args.strategy = args.strategy.lower()

    return args

def bc2int(bc):
    mapped = bc.translate(MAPPING).encode()
    bc_vec = np.frombuffer(mapped, dtype='uint8')
    return bc_vec

def compute_dist_matrix(barcodes):
    # barcode matrix is the base 5 representation of the barcodes.
    bc_len = max(map(len, barcodes))

    barcodes_matrix = np.array([bc2int(b) for b in barcodes])

    all_combinations = np.array([bc2int(''.join(bc)) for bc in product(NUCLS, repeat=bc_len)], dtype=np.uint8)

    distances = (cdist(barcodes_matrix, all_combinations, metric='hamming') * bc_len).astype(np.uint8)

    return distances.T

def main():
    '''
    - Open the 2 index files at the same time
    - Compare both fwd and reverse index with barcodes
    - If both match with less than [max_mismatches], it's a hit
        - If there are multiple hits, we take the best one.
        - If there are multiple best hits, we don't choose
    - We write (read_id, fwd_index, rev_index, sample_id) in a file
    '''

    args = parse_args()

    barcodes = pd.read_csv(args.meta, dtype=str,
                           names=["sample_name", "fwd_bc", "rev_bc"])
    barcodes.sample_name.str.replace("[^a-zA-Z0-9_.]", "_")

    distances = [compute_dist_matrix(barcodes.fwd_bc), compute_dist_matrix(barcodes.rev_bc)]

    if args.revMapping:
        distances = distances[::-1]

    mapping = barcodes.sample_name.values

    # Iterate over fastq and store results to speed up processing
    fwd_handle = FastqGeneralIterator(open(args.fwd, "r"))
    rev_handle = FastqGeneralIterator(open(args.rev, "r"))

    for (hdr_fwd, seq_fwd, _), (hdr_rev, seq_rev, _) in zip(fwd_handle, rev_handle):
        mismatches = sum(distances[i][int(seq.translate(MAPPING_BASE5), 5)]
                         for (i, seq) in enumerate([seq_fwd, seq_rev]))
        sample_idx = np.argwhere(mismatches <= 2*args.max_mismatches)[:, 0]

        # Default
        sample_names = ["_UNKNOWN_"]

        if len(sample_idx) == 1 or (len(sample_idx) > 1 and args.strategy == 'all'):
            sample_names = mapping[sample_idx]
        elif len(sample_idx) > 1:
            best_score = mismatches.min()
            best_matches = np.argwhere(mismatches == best_score)[:, 0]

            if args.strategy in ['best', 'min_all']:
                if len(best_matches) == 1:
                    sample_names = mapping[best_matches]
                else:
                    if args.strategy == 'min_all':
                        sample_names = mapping[best_matches]
            elif args.strategy == 'discard':
                pass
                
        id_fwd = hdr_fwd.split(" ")[0]
        id_rev = hdr_rev.split(" ")[0]

        if id_fwd == id_rev:
            for sample_name in sample_names:
                print("\t".join([id_fwd, seq_fwd, seq_rev, sample_name]))
        else:
            print("error: forward and reverse index don't match")
            exit(1)

    fwd_handle.close()
    rev_handle.close()

if __name__ == '__main__':
    main()
