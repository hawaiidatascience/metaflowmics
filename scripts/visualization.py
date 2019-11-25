import argparse

import pandas as pd
import numpy as np
from Bio import SeqIO

import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=str, default='.')
    parser.add_argument('--thresh', type=int, default=97)
    parser.add_argument('--show', action='store_true', default=False)    
    args = parser.parse_args()

    return args

class Loader:

    def __init__(self, args):
        self.root = args.root
        self.thresh = args.thresh
        self.path = {
            'fasta': f'{args.root}/sequences_{args.thresh}.fasta',
            'shared': f'{args.root}/abundance_table_{args.thresh}.shared',
            'tax': f'{args.root}/annotations_{args.thresh}.taxonomy',
            'db': f'{args.root}/abundance_table_{args.thresh}.database'
        }
        self.loaded = {'fasta': False, 'shared': False, 'tax': False, 'db': False}

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
            elif key == 'db':
                self.load_db()
                self.loaded['shared'] = True
                self.loaded['tax'] = True                
            else:
                print(f'Cannot load {key} (unknown format)')
                return
            self.loaded[key] = True
            
        
    def load_shared(self):
        shared = (pd.read_csv(self.path['shared'], sep='\t', dtype={'Group': str})
                  .set_index('Group')
                  .drop(['label', 'numOtus'], axis=1))
        shared.columns.name = 'OTU'        
        self.shared = shared

    def load_tax(self):
        tax = pd.read_csv(self.path['tax'], index_col=0, sep='\t', header=None, names=['OTU', 'Size', 'Taxonomy'])
        tax = tax.Taxonomy.str.split(';', expand=True)
        tax.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        tax = tax.replace('', np.nan).dropna(axis=1, how='all')
        tax.index.name = 'OTU'

        self.tax = tax

    def load_db(self):
        db = (pd.read_csv(self.path['db'], index_col=0, sep='\t',
                          usecols=lambda x: x not in ['repSeq', 'repSeqName'])
              .rename(columns={'OTUConTaxonomy': 'tax'}))
        db.index.name = 'OTU'
        self.shared = db.drop('tax', axis=1).astype(int).T

        ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
        self.tax = (db.tax
                    .str.strip(';').str.split(';', expand=True)
                    .rename(columns={i: rank for (i, rank) in enumerate(ranks)}))
        
    def load_fasta(self):
        sequences = {seq.id: str(seq.seq) for seq in SeqIO.parse(self.path['fasta'], 'fasta')}
        
        self.fasta = sequences

    def set_prevalence(self, min_abund=0):
        self.load('db')
        self.prevalence = (self.shared > min_abund).mean(axis=0)

    def set_abundance(self):
        self.load('shared')
        self.abundance = self.shared.sum(axis=0)
        

def phylum_scatter(loader, rank='Phylum', select=None, min_group_size=5, show=False):
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

def main():
    '''
    '''

    args = parse_args()
    data_loader = Loader(args)

    # data_loader.load('shared')
    # data_loader.load('tax')
    data_loader.load('db')    
    
    data_loader.set_prevalence()
    data_loader.set_abundance()
    
    biclustering(data_loader, show=args.show)
    phylum_scatter(data_loader, show=args.show)
    phylum_scatter(data_loader, rank='Class', select=('Phylum', 'Proteobacteria'), show=args.show)

if __name__ == '__main__':
    main()
