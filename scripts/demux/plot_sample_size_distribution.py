'''
Required packages: glob, matplotlib

- Utility to help choosing the subsampling threshold
for 16S (potentially other markers as well) data analysis

- Extracts the number of reads in each samples and 
plot the distribution of the log counts.

'''

from math import log10
from glob import glob

import numpy as np
import pandas as pd

from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt

def count_samples(folder='.'):
    '''
    Count the number of reads in all forward fastq files in [folder]
    '''

    files = glob("{}/sample_counts*.csv".format(folder))
    summaries = pd.concat([pd.read_csv(summary, index_col=0, dtype=str, header=None)
                           for summary in files],
                          axis=1, sort=True)

    return summaries.fillna(0).astype(int).sum(axis=1).values

def plot(counts, display=False):
    '''
    Display histogram of log counts
    '''

    fig, ax = plt.subplots()
    ax.hist(np.log10(counts), bins=max(5, len(counts)//10))
    ax.set_xlabel('Sample size', fontsize=14)
    ax.set_ylabel('Occurrence', fontsize=14)
    ax.xaxis.set_major_locator(MultipleLocator(1))

    minor_ticks = [log10(j*10**i) for j in range(2,10) for i in ax.get_xticks()]
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(["{:,}".format(int(10**i)) for i in ax.get_xticks()])

    if not display:
        fig.savefig('sample_sizes.pdf', transparent=True)
    else:
        plt.show()

def main():
    '''
    Argument: folder where the demultiplexed samples are located
    The reads need to be formatted as .fastq (not gzipped)
    '''

    counts = count_samples()
    plot(counts)

if __name__ == '__main__':
    main()
