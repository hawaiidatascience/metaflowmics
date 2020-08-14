import re
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.models.ranges import FactorRange
from bokeh.layouts import gridplot
from bokeh.palettes import Turbo256, linear_palette
from bokeh.models import ColorBar, LinearColorMapper, BasicTicker, Legend
from bokeh.transform import transform

RANKS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
TOOLS = ['hover', 'box_zoom', 'reset']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--shared', type=str)
    parser.add_argument('-t', '--taxonomy', type=str)
    parser.add_argument('--thresh', type=str)
    parser.add_argument('--min-abund', type=float, default=0.01)
    parser.add_argument('--max-samples', type=int, default=500)        

    args = parser.parse_args()

    return args

def load_shared(path):
    header = next(open(path)).split('\t')
    dtypes = dict((col, str) if col in ['Group'] else (col, int)
                  for col in header)
    data = pd.read_csv(path, sep='\t', dtype=dtypes,
                       keep_default_na=False, low_memory=False,
                       usecols=lambda x: x not in ['label', 'numOtus'])
    return data.set_index('Group')

def load_tax(path):
    df = (pd.read_csv(path, index_col='OTU', sep='\t').drop('Size', axis=1).Taxonomy
          .str.strip(';')
          .str.replace('-', '_', regex=False)
          .str.replace('[\(\)]', '')
          .str.split(';', expand=True))
    df.columns = RANKS[:df.shape[1]]

    return df

def data(args):
    return {'tax': load_tax(args.taxonomy),
            'shared': load_shared(args.shared)}

def group_otus(shared, tax, rank, min_abund=0.01):
    data = shared.T.groupby(tax[rank]).agg('sum').rename_axis(columns='sample', index=rank)

    stack_data = data.stack().rename('proportion').reset_index()
    lim = (data.sum()*min_abund).loc[stack_data['sample']].to_numpy()

    filler = '_Others (< {:.0%})'.format(min_abund)

    stack_data.loc[stack_data.proportion < lim, rank] = filler

    data = stack_data.groupby(['sample', rank]).proportion.agg(sum).unstack().fillna(0)

    return data

def scatter(shared, tax, rank='Phylum', output='scatter.html'):
    if rank == 'Class':
        subset = tax.Phylum == 'Proteobacteria'
    else:
        subset = ~tax.Phylum.isnull()

    data = pd.DataFrame({
        'abundance': shared.loc[:, subset].sum(),
        'prevalence': (shared.loc[:, subset] > 0).sum(),
    })
    data = pd.concat([data, tax.loc[data.index]], axis=1)
    data.index.name = 'OTU'

    otu_groups = tax.loc[data.index].reset_index().groupby(rank).OTU.agg(list)

    tooltips = list(zip(tax.columns, '@'+tax.columns)) + [('abundance', '@abundance{0,0}'), ('prevalence', '@prevalence')]

    output_file(output)

    plots = []
    for name, otus in otu_groups.iteritems():
        p = figure(title=name.capitalize(), tooltips=tooltips, tools=TOOLS)
        p.circle(x='prevalence', y='abundance', size=5, alpha=0.5, hover_color='red', source=data.loc[otus])
        p.xaxis.axis_label = 'OTU prevalence'
        p.yaxis.axis_label = 'OTU abundance'
        p.title.text_font_size = '14pt'
        plots.append(p)

    grid = gridplot(plots, ncols=5, plot_width=300, plot_height=300)
    save(grid)


def stackplot(shared, tax, rank='Phylum', output='barplot.html', min_abund=0.01, offset=50):

    x_sums = shared.sum(axis=1)
    tooltips = [('sample', '@sample'),
                (rank, '@'+rank),
                ('proportion', '@proportion{0.0%}'),
                ('total abundance', '@total_abundance')]

    if rank == 'Class':
        shared = shared.loc[:, tax.Phylum == 'Proteobacteria']
        tax = tax.loc[shared.columns]
        tooltips += [('total_proteobacteria', '@total_proteobacteria')]
        prot_sums = shared.sum(axis=1)

    # Normalize to ratio to better visualize distribution
    shared = (shared.T / shared.sum(axis=1)).T

    # Group by rank
    shared = group_otus(shared, tax, rank, min_abund=min_abund)

    hue_order = shared.sum().sort_values(ascending=False).index
    
    filler = '_Others (< {:.0%})'.format(min_abund)
    if filler in hue_order:
        hue_order = hue_order.drop(filler).append(pd.Index([filler]))
    cmap = dict(zip(hue_order, linear_palette(Turbo256, len(hue_order))))
    cmap[filler] = '#cccccc'
    
    # If it's the first plot to overlay
    p = figure(y_range=FactorRange(*sorted(shared.index, reverse=True)),
               plot_height=2*offset+shared.shape[0]*15,
               plot_width=2*offset+1000,
               title=Path(output).stem.replace('_', ' '),               
               min_border=offset,
               y_axis_label='Sample', x_axis_label='Taxonomic composition',
               tooltips=tooltips)

    bar_data = pd.DataFrame({
        'left': shared[hue_order].shift(axis=1).cumsum(axis=1).fillna(0).stack(),
        'right': shared[hue_order].cumsum(axis=1).stack(),
        'proportion': shared[hue_order].stack()
    }).swaplevel(0, 1)

    samples = bar_data.index.get_level_values('sample')

    bar_data['color'] = [cmap[otu] for otu in bar_data.index.get_level_values(rank)]
    bar_data['total_abundance'] = x_sums.map(lambda x: f'{x:,}').loc[samples].to_numpy()

    if rank.lower() == 'class':
        bar_data['total_proteobacteria'] = prot_sums.map(lambda x: f'{x:,}').loc[samples].to_numpy()

    # plot each level one by one
    for i, level in enumerate(hue_order):
        data_i = bar_data.loc[level].assign(**{rank: level})

        p.hbar(left='left', right='right', height=0.8,
               y='sample', color='color',
               line_color='black', line_width=1.2,
               source=data_i.reset_index(),
               legend_label=level,
               name=level)

    leg_items = p.legend.items.copy()
    p.legend.items.clear()
    
    legend = Legend(items=leg_items,
                    padding=0, spacing=0, border_line_color=None)

    p.add_layout(legend, 'right')
        
    output_file(output)
    save(p)
    
def preproc_for_heatmap(shared, tax, rank, min_abund=0.01):

    data = group_otus(shared, tax, rank, min_abund=min_abund)

    # z-score first and biclustering
    z_data = (data - data.mean()) / data.std()

    # hierarchical clustering on both rows and columns
    try:
        sample_links = leaves_list(linkage(z_data, method='average', metric='braycurtis',
                                        optimal_ordering=True))
    except ValueError:
        print('Something went wrong with the hierarchical clustering on samples')
        sample_links = np.arange(z_data.shape[0])
    try:
        feature_links = leaves_list(linkage(z_data.T, method='average', metric='braycurtis',
                                        optimal_ordering=True))
    except ValueError:
        print('Something went wrong with the hierarchical clustering on features')
        feature_links = np.arange(z_data.shape[1])

    clustered_index = pd.MultiIndex.from_product([data.index[sample_links],
                                                  data.columns[feature_links]])
        
    data = pd.DataFrame({'z_score': z_data.stack(), 'abundance': data.stack()}).loc[clustered_index]

    tax[rank] = tax[rank].map(lambda x: float('nan')
                              if re.search('uncultured|unclassified|incertae', x)
                              else x)
    tax_cols = tax.columns[:tax.columns.get_loc(rank)]

    tax = tax.dropna(how='any').groupby(rank, as_index=False)[tax_cols].agg('first')

    data = data.reset_index().merge(tax, how='left', left_on=rank, right_on=rank)

    return data

def clustermap(data, fig_dir='.', output='heatmap.html'):

    tooltips = dict(zip(data.columns, '@' + data.columns))
    (samples, features) = [data[col] for col in data.columns[:2]]

    p = figure(plot_height=10*len(samples.unique()),
               plot_width=10*len(features.unique()),
               title=Path(output).stem.replace('_', ' '),
               x_range=features.unique(),
               y_range=samples.unique().tolist(),
               x_axis_location='above',
               tools=TOOLS, tooltips=list(tooltips.items()))

    mapper = LinearColorMapper(palette='Viridis256',
                               low=data.z_score.min(),
                               high=data.z_score.max())

    p.rect(x=features.name, y=samples.name, width=1, height=1, source=data.reset_index(),
           line_color=None, fill_color=transform('z_score', mapper))

    color_bar = ColorBar(color_mapper=mapper, formatter=p.xaxis.formatter,
                         ticker=BasicTicker(desired_num_ticks=5),
                         location=(0,0))
    
    p.add_layout(color_bar, 'right')

    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = '8pt'
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = 'vertical'
    p.title.text_font_size = '14pt'

    output_file(output)
    save(p)
    
if __name__ == '__main__':
    args = parse_args()

    data = data(args)

    scatter(rank='Phylum', output='scatter_OTU{}_by-phylum.html'.format(args.thresh),
            **data.copy())
    scatter(rank='Class', output='scatter_OTU{}_proteobacteria_by-class.html'.format(args.thresh),
            **data.copy())

    suffix = ''
    if data['shared'].shape[0] > args.max_samples:
        data['shared'] = data['shared'].sample(args.max_samples).sort_index()
        suffix = '_minAbund-{}_{}-random-samples'.format(args.min_abund, args.max_samples)

    stackplot(rank='Phylum', min_abund=args.min_abund,
              output='barplot_OTU{}_by-Phylum{}.html'.format(args.thresh, suffix),
              **data.copy())
    stackplot(rank='Class', min_abund=args.min_abund,
              output='barplot_OTU{}_preoteobacteria_by-class{}.html'.format(args.thresh, suffix),
              **data.copy())

    cluster_data = preproc_for_heatmap(**data, rank='Order', min_abund=0.05)
    clustermap(cluster_data, output='biclutered_heatmap_sample-by-order-{}{}.html'
               .format(args.thresh, suffix))
    
    cluster_data = preproc_for_heatmap(**data, rank='Family', min_abund=0.05)
    clustermap(cluster_data, output='biclutered_heatmap_sample-by-family-{}{}.html'
               .format(args.thresh, suffix))
