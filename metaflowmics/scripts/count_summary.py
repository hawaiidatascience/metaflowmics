import argparse
from pathlib import Path

import pandas as pd

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--list', type=str)
    parser.add_argument('--count', type=str)    
    args = parser.parse_args()

    return args

def count(list_file, count_file):
    ct = pd.read_csv(count_file, index_col=0, sep='\t',
                     usecols=lambda x: x != 'total')

    otus = pd.read_csv(list_file, sep='\t', usecols=lambda x: x not in ['label', 'numOtus']).T
    otus = otus[0].str.split(',').explode().reset_index().set_index(0)['index']

    ct.index = otus.loc[ct.index]
    ct = ct.groupby(level=0).agg(sum)
    
    nseqs = ct.sum(axis=0)
    sobs = (ct > 0).sum(axis=0)

    summary_file = pd.DataFrame({'nseqs': nseqs, 'sobs': sobs})
    summary_file.index.name = 'group'

    summary_file.to_csv('{}.groups.summary'.format(Path(list_file).stem), sep='\t')

def main():
    '''
    '''

    args = parse_args()
    count(args.list, args.count)

if __name__ == '__main__':
    main()
