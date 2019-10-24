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
    parser.add_argument('--fwd', type=str)
    parser.add_argument('--rev', type=str)
    parser.add_argument('--mapping', type=str)
    args = parser.parse_args()

    return args

def main():
    '''
    - Read demux_info file, and both read files in a synchronized way
    - Write read files to a file depending on the sample assignment in demux_info
    '''

    args = parse_args()

    demux_info = pd.read_csv(args.mapping, header=None, index_col=0, sep="\t", dtype=str,
                             names=['bc_fwd', 'bc_rev', 'sample_name'])

    print('Preparing handles.')
    handles = {}
    for sample in demux_info['sample_name'].unique():
        if sample != '_UNKNOWN_':
            handles[sample+'fwd'] = open('{}_R1.fastq'.format(sample), 'w')
            handles[sample+'rev'] = open('{}_R2.fastq'.format(sample), 'w')

    parser_fwd = FastqGeneralIterator(open(args.fwd, 'r'))
    parser_rev = FastqGeneralIterator(open(args.rev, 'r'))

    print('Starting demultiplexing')
    for seq_nb, (seq_fwd, seq_rev) in enumerate(zip(parser_fwd, parser_rev)):
        id_fwd = seq_fwd[0].split(' ')[0]
        id_rev = seq_rev[0].split(' ')[0]

        if id_fwd != id_rev:
            print("Sequence #{}: {} and {} do not match"
                  .format(seq_nb, id_fwd, id_rev))
            continue

        sample_assignment = demux_info.loc[id_fwd, "sample_name"]

        if sample_assignment == '_UNKNOWN_':
            continue

        handles[sample_assignment+'fwd'].write('@{}\n{}\n+\n{}\n'.format(*seq_fwd))
        handles[sample_assignment+'rev'].write('@{}\n{}\n+\n{}\n'.format(*seq_rev))

    for sample in demux_info['sample_name'].unique():
        if sample != '_UNKNOWN_':
            handles[sample+'fwd'].close()
            handles[sample+'rev'].close()

    print("Demultiplexing finished.")

if __name__ == '__main__':
    main()
