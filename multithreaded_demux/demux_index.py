#!/usr/bin/env python3

'''
- Read index and barcodes
- Match each read_id to a sample_id
- Write the result in a file
'''

import argparse

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np

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
    parser.add_argument('--strategy', type=str, help='Strategy to handle multimappers')    
    args = parser.parse_args()

    args.strategy = args.strategy.lower()

    return args

def get_mismatches(query, targets):
    '''
    Compare a query sequence with all barcodes in targets
    '''
    mismatches = np.zeros(len(targets), dtype=np.uint8)

    for i, char in enumerate(query):
        mismatches[targets[:, i] != char] += 1

    return mismatches

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

    # Reversed forward and reverse index
    bc_seq = [
        np.array(barcodes.fwd_bc.apply(list).tolist()),
        np.array(barcodes.rev_bc.apply(list).tolist())
    ]

    if args.revMapping:
        bc_seq = bc_seq[::-1]

    mapping = barcodes.sample_name.values

    # Iterate over fastq and store results to speed up processing
    fwd_handle = FastqGeneralIterator(open(args.fwd, "r"))
    rev_handle = FastqGeneralIterator(open(args.rev, "r"))

    seen = [{}, {}]

    for (hdr_fwd, seq_fwd, _), (hdr_rev, seq_rev, _) in zip(fwd_handle, rev_handle):
        both_mismatches = []

        for i, seq in enumerate([seq_fwd, seq_rev]):
            # Check if we already saw this index
            mismatches = seen[i].get(seq, None)

            if mismatches is None:
                mismatches = get_mismatches(seq, bc_seq[i])
                seen[i][seq] = mismatches

            both_mismatches.append(mismatches)

        mismatches = sum(both_mismatches)
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
            else:
                print('ERROR: Unknown strategy {}. Aborting.'.format(args.strategy))
                exit(1)
                
        id_fwd = hdr_fwd.split(" ")[0]
        id_rev = hdr_rev.split(" ")[0]

        if id_fwd == id_rev:
            for sample_name in sample_names:
                print("\t".join([id_fwd, seq_fwd, seq_rev, sample_name]), end="\n\r")
        else:
            print("error: forward and reverse index don't match")
            exit(1)

    fwd_handle.close()
    rev_handle.close()

if __name__ == '__main__':
    main()
