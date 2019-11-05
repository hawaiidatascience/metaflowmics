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
    parser.add_argument('--meta', type=str, help='Path to metadata')
    parser.add_argument('--max_mismatches', type=int, default=1)
    parser.add_argument('--strategy', type=str, choices=['best', 'min_all', 'discard'], default='best', help='Strategy to handle multimappers')    
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

    barcodes = pd.read_csv(args.meta, dtype=str, names=["sample_name", "fwd_bc"])
    barcodes.sample_name.str.replace("[^a-zA-Z0-9_.]", "_")

    # Reversed forward and reverse index
    bc_seq = np.array(barcodes.fwd_bc.apply(list).tolist())

    mapping = barcodes.sample_name.values

    # Iterate over fastq and store results to speed up processing
    handle = FastqGeneralIterator(open(args.fwd, "r"))

    seen = {}

    for (hdr, seq, _) in handle:
        # Check if we already saw this index
        mismatches = seen.get(seq, None)

        if mismatches is None:
            mismatches = get_mismatches(seq, bc_seq)
            seen[seq] = mismatches

        sample_idx = np.argwhere(mismatches <= args.max_mismatches)[:, 0]

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
                
        idx = hdr.split(" ")[0]

        for sample_name in sample_names:
            print("\t".join([idx, seq, '', sample_name]))

    handle.close()

if __name__ == '__main__':
    main()
