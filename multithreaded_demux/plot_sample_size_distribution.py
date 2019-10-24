'''
Required packages: glob, matplotlib

- Utility to help choosing the subsampling threshold
for 16S (potentially other markers as well) data analysis

- Extracts the number of reads in each samples and 
plot the distribution of the log counts.

'''

import sys
from itertools import chain
from math import log10
import subprocess as sp
from glob import glob

from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt

def count_samples(folder='.'):
    '''
    Count the number of reads in all forward fastq files in [folder]
    '''
    
    nb_reads = []
    for f in glob("{}/*_R1*.fastq".format(folder)):
        parsed_sys_call = sp.check_output(['wc', '-l', f]).decode().split(' ')[0]
        log_count = log10(int(parsed_sys_call))
        nb_reads.append(log_count)
    
    return nb_reads

def plot(counts, display=False):
    '''
    Display histogram of log counts
    '''

    fig, ax = plt.subplots()
    ax.hist(counts, bins=30)
    ax.set_xlabel('Sample size', fontsize=14)
    ax.set_ylabel('Occurrences', fontsize=14)
    ax.xaxis.set_major_locator(MultipleLocator(1))

    minor_ticks = list(chain(*[[log10(j*10**i) for j in range(2,10)] for i in ax.get_xticks()]))
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(["{:,}".format(int(10**i)) for i in ax.get_xticks()])

    if not display:
        fig.savefig('sample_sizes.png')
    else:
        plt.show()

def main():
    '''
    Argument: folder where the demultiplexed samples are located
    The reads need to be formatted as .fastq (not gzipped)
    '''

    folder = sys.argv[1]
    counts = count_samples(folder)
    plot(counts)

if __name__ == '__main__':
    main()
