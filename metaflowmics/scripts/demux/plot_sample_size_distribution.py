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

from bokeh.plotting import figure
from bokeh.io import save, output_file
from bokeh.models import tickers

from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt

def count_samples(folder='.'):
    '''
    Count the number of read pair/sample in [folder]
    '''

    files = glob("{}/sample_counts*.csv".format(folder))
    summaries = pd.concat([pd.read_csv(summary, index_col=0, dtype=str, header=None)
                           for summary in files],
                          axis=1, sort=True)

    return summaries.fillna(0).astype(int).sum(axis=1)

def describe_bokeh_obj(x):
    info = ['{:50}: {}'.format(x,y) for (x,y) in x.properties_with_values().items()]

    print('\n'.join(info))

def plot_bokeh(counts):

    hist, edges = np.histogram(np.log10(counts), bins=max(5, len(counts)//10), density=False)
    hist_df = pd.DataFrame({'count': hist,
                            "left": edges[:-1],
                            "right": edges[1:]})
    hist_df["interval"] = ["{:,} - {:,}".format(int(10**left), int(10**right))
                           for left, right in zip(hist_df["left"], hist_df["right"])]

    x_min = int(min(edges))
    x_max = max(4, 1+int(max(edges)))
    
    p = figure(plot_height=800, plot_width=800,
               x_range=[x_min, x_max], tools='hover,box_zoom',
               tooltips=[('Size range', '@interval'),
                         ('#Samples in interval', str("@count"))],
                title="Sample size distribution",
                x_axis_label="Sample read count",
                y_axis_label="Occurrences")

    p.quad(bottom=0, top="count", left="left", 
           right="right", source=hist_df, fill_color="SteelBlue", 
           line_color="black", fill_alpha=0.7,
           hover_fill_alpha=1.0, hover_fill_color="Tan")

    ticks = list(range(x_min, x_max))
    minor_ticks = np.log10([i*10**j for i in range(1, 10) for j in ticks])
    
    p.xaxis.ticker = tickers.FixedTicker(ticks=ticks, minor_ticks=minor_ticks)
    p.xaxis.major_label_overrides = {tick: "{:,}".format(int(10**tick)) for tick in ticks}
    p.yaxis.minor_tick_line_color = None

    p.axis.major_label_text_font_size = "12pt"
    p.axis.axis_label_text_font_size = "14pt"
    p.title.text_font_size = '18pt'

    output_file('sample_sizes.html')
    save(p)

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
    counts.to_csv('sample_sizes.csv', header=False)
    plot_bokeh(counts)

if __name__ == '__main__':
    main()
