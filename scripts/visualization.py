import argparse
from collections import defaultdict

import pandas as pd
import numpy as np
from Bio import SeqIO

import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")

def parse_args():
    '''
    Retrieve all possible arguments from the command line
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=str, default='.', descr='Directory with all pipeline outputs (.shared, .taxonomy or .database)')
    parser.add_argument('--meta', type=str, default='', descr='Path to metadata (optional)')
    parser.add_argument('--thresh', type=int, default=100, descr='OTU threshold')
    parser.add_argument('--show', action='store_true', default=False, descr='Display the figure. If not set, saves it as PDF')
    parser.add_argument('--skip-meta', action='store_true', default=False, descr='Skip plots including metadata information')
    args = parser.parse_args()

    return args

class Loader:
    '''
    Data loader object to facilitate loading and organizing data
    '''

    def __init__(self, args):
        self.root = args.root
        self.thresh = args.thresh
        self.path = {
            'fasta': f'{args.root}/sequences_{args.thresh}.fasta',
            'shared': f'{args.root}/abundance_table_{args.thresh}.shared',
            'tax': f'{args.root}/annotations_{args.thresh}.taxonomy',
            'db': f'{args.root}/abundance_table_{args.thresh}.database',
            'meta': args.meta
        }
        self.loaded = defaultdict(bool)

    def load(self, key):
        if self.loaded[key]:
            return
        else:
            print('Loading {}'.format(key))
            if key == 'fasta':
                self.load_fasta()
            elif key == 'shared':
                self.load_shared()
            elif key == 'tax':
                self.load_tax()
            elif key == 'meta':
                self.load_meta()
            elif key == 'db':
                self.load_db()
                self.loaded['shared'] = True
                self.loaded['tax'] = True                
            else:
                print(f'Cannot load {key} (unknown format)')
                return
            self.loaded[key] = True
        
    def load_shared(self):
        '''
        Loader for mothur shared file
        '''
        
        shared = (pd.read_csv(self.path['shared'], sep='\t', dtype={'Group': str})
                  .set_index('Group')
                  .drop(['label', 'numOtus'], axis=1))
        shared.columns.name = 'OTU'
        shared.index.name = 'SampleID'
        self.shared = shared

    def load_tax(self):
        '''
        Loader for mothur taxonomy file
        '''

        tax = pd.read_csv(self.path['tax'], index_col=0, sep='\t', header=None, names=['OTU', 'Size', 'Taxonomy'])
        tax = tax.Taxonomy.str.split(';', expand=True)
        tax.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        tax = tax.replace('', np.nan).dropna(axis=1, how='all')
        tax.index.name = 'OTU'

        self.tax = tax

    def load_db(self):
        '''
        Loader for mothur database file
        '''

        db = (pd.read_csv(self.path['db'], index_col=0, sep='\t',
                          usecols=lambda x: x not in ['repSeq', 'repSeqName'])
              .rename(columns={'OTUConTaxonomy': 'tax'}))
        db.index.name = 'OTU'

        self.shared = db.drop('tax', axis=1).astype(int).T
        self.shared.index.name = 'SampleID'        

        ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
        self.tax = (db.tax
                    .str.strip(';').str.split(';', expand=True)
                    .rename(columns={i: rank for (i, rank) in enumerate(ranks)}))

    def load_meta(self, path=None):
        '''
        Loader for metadata file
        '''
        
        if not self.path['meta']:
            return

        self.meta = pd.read_csv(self.path['meta'], dtype=str)
        
    def load_fasta(self):
        '''
        Loader for fasta file
        '''
        
        sequences = {seq.id: str(seq.seq) for seq in SeqIO.parse(self.path['fasta'], 'fasta')}
        
        self.fasta = sequences

    def set_prevalence(self, min_abund=0):
        '''
        Computes prevalence
        '''
        
        self.load('db')
        self.prevalence = (self.shared > min_abund).mean(axis=0)

    def set_abundance(self):
        '''
        Computes abundance of each OTU
        '''
        
        self.load('shared')
        self.abundance = self.shared.sum(axis=0)
        

def phylum_scatter(loader, rank='Phylum', select=None, min_group_size=5, show=False):
    '''
    Abundance vs Prevalence scatter plot for each OTU. 
    Each facet corresponds to one of the possible ranks at the given level (default is Phylum)
    '''
    
    summaries = pd.DataFrame(
        {'Prevalence': loader.prevalence,
         'Abundance': loader.abundance,
         rank: loader.tax.loc[loader.abundance.index, rank]})

    if select is not None:
        rank_s, label = select
        summaries = summaries[loader.tax[rank_s] == label]

    # Filter out Phyla that are not well represented
    others = summaries.reset_index().groupby(rank).filter(lambda x: len(set(x['OTU'])) <= min_group_size)[rank].unique()
    summaries.loc[np.isin(summaries[rank], others), rank] = '_Other(<={} OTUs)'.format(min_group_size)

    col_wrap = int(2+np.sqrt(len(summaries[rank].unique())))
    
    sns.set(font_scale=0.8)
    g = sns.FacetGrid(data=summaries, col=rank, col_wrap=col_wrap, sharex=False, hue=rank,
                      palette='Set1', col_order=sorted(summaries[rank].unique()))
    g.map(plt.scatter, 'Abundance', 'Prevalence', alpha=0.5, s=40)
    g.despine(left=True)
    g.set_titles("{col_name}", fontweight='bold', fontsize=14)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.05)

    plt.savefig("scatterplot_{}_abd-prev_{}.pdf".format(rank, loader.thresh), transparent=True)

    if show:
        plt.show()

def biclustering(loader, top=100, cols=['Phylum', 'Class', 'Order'], show=False):
    '''
    Biclustering of sample x OTU abundance matrix (z-scores). Display {cols} on the right side.
    '''

    sorted_otus = (loader.abundance * loader.prevalence).sort_values(ascending=False).index
    selected_otus = sorted_otus[loader.shared.std() > 0][:top]

    matrix = loader.shared.T.loc[selected_otus]
    taxa = loader.tax.loc[selected_otus, cols]
    slen = taxa.applymap(len).max()

    for i, col in enumerate(cols):
        taxa[col] = taxa[col].str.pad(slen[col]+2, side='right', fillchar=' ')
        formatter = '{:' + str(slen[col]) + '}'
        cols[i] = formatter.format(col)

    matrix.index = np.sum(taxa, axis=1)

    sns.set(font_scale=0.5)
    g = sns.clustermap(matrix, figsize=(20, 15), row_colors=None, z_score=0)
    plt.subplots_adjust(bottom=0.2, right=0.7)

    ax = g.ax_heatmap
    right_clabel_pos = ax.get_xmajorticklabels()[-1]._x

    ax.text(5 + right_clabel_pos, 0, '  '.join(cols), fontsize=5, fontweight='bold')

    plt.savefig("biclustering_{}.pdf".format(loader.thresh), transparent=True)

    if show:
        plt.show()

def stacked_bars(loader, level='Phylum', factors=None, norm=True, n_top=-1, out_prefix='bars'):
    '''
    Stacked bar graph showing the levels (default is Phylum).
    - x-axis is factors[0] (or sample id if empty)
    - row facets are factors[1]
    - col facets are factors[2]
    If there are too many possible annotations, n_top can limit the annotations to the {n_top} most abundance
    '''
    
    shared = loader.shared.T.copy()
    shared[level] = loader.tax[level]
    shared_by_tax = shared.groupby(level).agg(sum).T

    if n_top > 0: 
        topn_list = shared_by_tax.sum().sort_values(ascending=False).index
        shared_by_tax = shared_by_tax[topn_list]

    if norm:
        shared_by_tax = shared_by_tax / shared_by_tax.sum(axis=1)[:, None]
        out_prefix += '_norm'

    if factors is not None:
        meta = loader.meta[factors].reindex(shared_by_tax.index)

        shared_by_tax = pd.concat([shared_by_tax, meta], sort=True, axis=1)
        shared_by_tax = shared_by_tax.dropna().groupby(factors).agg(np.mean)

    else:
        factors = ['SampleID']

    if shared_by_tax.size == 0:
        print('No samples with all metadata. Aborting')
        return

    plot_prms = {}
    if len(factors) >= 2:
        plot_prms['row'] = factors[1]
    if len(factors) == 3:
        plot_prms['col'] = factors[2]

    all_colors = sns.hls_palette(shared_by_tax.shape[1], l=.4, s=.7)
    palette = {tax: color for color, tax in zip(all_colors, shared_by_tax.columns)}

    shared_by_tax = shared_by_tax.reset_index().melt(id_vars=factors)

    g = sns.FacetGrid(data=shared_by_tax, sharey=False, **plot_prms)
    g.map(stacked_single, factors[0], 'value', level,
          palette=palette, x_values=sorted(shared_by_tax[factors[0]].unique()))

    plt.savefig(f'{out_prefix}_stacked_bars.pdf', transparent=True)    

def stacked_single(x, y, hue, palette=None, x_values=None, **kwargs):
    '''
    Sub-routine to plot stacked bar graph
    '''

    ax = plt.gca()
    df = pd.DataFrame({'x': x, 'y': y, 'hue': hue}).pivot('x', 'hue')
    df = df.reindex(index=x_values).fillna(0)
    df.columns = df.columns.get_level_values('hue')

    prev_bars = np.zeros(df.shape[0])
    width = 5 / len(x_values)

    for (label, bar) in df.T.iterrows():
        ax.bar(bar.index, bar.values, color=palette[label],
               edgecolor='k', width=width, linewidth=0.2,
               bottom=prev_bars, label=label)
        prev_bars += bar.values
    
    
def main():
    '''
    Script runner
    '''

    args = parse_args()
    data_loader = Loader(args)

    data_loader.load('db')
    data_loader.load('meta')
    
    data_loader.set_prevalence()
    data_loader.set_abundance()

    biclustering(data_loader, show=args.show)
    phylum_scatter(data_loader, show=args.show)
    phylum_scatter(data_loader, rank='Class', select=('Phylum', 'Proteobacteria'), show=args.show)

    if not args.skip_metadata:
        stacked_bars(data_loader, level='Phylum', norm=True, n_top=-1, out_prefix='bars')
    
if __name__ == '__main__':
    main()
