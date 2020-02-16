#!/usr/bin/env python3

'''
Use the demux_info file from demux_index
and sort reads according to their sample assignments.

For large number of samples, you might need to increased the 
allow number of opened files (e.g. `ulimit -n 4096`)
'''

import argparse

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd

def parse_args():
    '''
    Get command line arguments
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--fastqs', type=str, nargs='+')
    parser.add_argument('--mapping', type=str)
    args = parser.parse_args()

    return args

def main():
    '''
    - Read demux_info file, and both read files in a synchronized way
    - Write read files to a file depending on the sample assignment in demux_info
    '''

    args = parse_args()

    demux_info = pd.read_csv(args.mapping, header=None, index_col=0, sep="\t", dtype=str).dropna(axis=1, how='all')
    index_orient = ['fwd', 'rev'][:demux_info.shape[1]-3]
    demux_info.columns = ['rid'] + index_orient + ['sample_name', 'mismatches']

    read_orient = ['fwd', 'rev'][:len(args.fastqs)]
    
    print('Preparing handles.')
    handles = {}
    for sample in demux_info['sample_name'].unique():
        if not pd.isnull(sample):
            for i, orient in enumerate(read_orient, 1):
                handles[sample+orient] = open('{}_R{}.fastq'.format(sample, i), 'w')

    parsers = [FastqGeneralIterator(open(fastq, 'r')) for fastq in args.fastqs]

    print('Starting demultiplexing')
    for seq_nb, sequences in enumerate(zip(*parsers)):
        ids = [seq[0].split()[0] for seq in sequences]

        if len(ids) > 1:
            if ids[0] != ids[1]:
                print("Sequence #{}: {} (fwd) and {} (rev) do not match. The forward and reverse read files seem to be out of order"
                      .format(seq_nb, *ids))
                exit(42)

        sample_assignment = demux_info.loc[ids[0], "sample_name"]

        if pd.isnull(sample_assignment):
            continue

        for orient, seq in zip(read_orient, sequences):
            handles[sample_assignment + orient].write('@{}\n{}\n+\n{}\n'.format(*seq))

    for sample in demux_info['sample_name'].unique():
        if not pd.isnull(sample):
            for orient in read_orient:
                handles[sample + orient].close()

    print("Demultiplexing finished.")

if __name__ == '__main__':
    main()
