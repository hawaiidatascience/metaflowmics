import argparse

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.palettes import Dark2, Set1, Set2, Set3, Category20
from bokeh.models import ColorBar, LinearColorMapper, BasicTicker
from bokeh.transform import transform

RANKS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
TOOLS = ['hover', 'box_zoom', 'reset']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--shared', type=str)
    parser.add_argument('-t', '--taxonomy', type=str)
    parser.add_argument('--thresh', type=str)
    parser.add_argument('--ntop', type=int, default=100)    
    
    args = parser.parse_args()

    return args

def load_shared(path):
    df = pd.read_csv(path, index_col='Group', sep='\t').drop(['label', 'numOtus'], axis=1)
    return df

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

def group_otus(shared, tax, rank, top=None):
    data = shared.T.groupby(tax[rank]).sum().T

    suffix = ""
    if top is not None:
        if data.shape[1] > top:
            suffix = "_top-{}".format(top)
        top_otus = data.sum().sort_values(ascending=False).index[:top]
        data = data.loc[:, top_otus]

    return (data, suffix)

def scatter(shared, tax, thresh=100, rank='Phylum'):
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
    
    output_file("scatter_OTU{}_by-{}.html".format(thresh, rank))

    plots = []
    for name, otus in otu_groups.iteritems():
        p = figure(title=name.capitalize(), tooltips=tooltips, tools=TOOLS)
        p.circle(x="prevalence", y="abundance", size=5, alpha=0.5, hover_color="red", source=data.loc[otus])
        p.xaxis.axis_label = 'OTU prevalence'
        p.yaxis.axis_label = 'OTU abundance'
        p.title.text_font_size = '14pt'
        plots.append(p)

    grid = gridplot(plots, ncols=5, plot_width=300, plot_height=300)
    save(grid)

def stackplot(shared, tax, thresh=100, rank='Phylum', top=20):

    if rank == 'Class':
        shared_sub = shared.loc[:, tax.Phylum == 'Proteobacteria']
        tax_sub = tax.loc[shared_sub.columns]
    else:
        shared_sub = shared.copy()
        tax_sub = tax.copy()

    (shared_sub, suffix) = group_otus(shared_sub, tax_sub, rank, top=top)
    outname = "barplot_OTU-{}_by-{}{}.html".format(thresh, rank, suffix)

    shared_sub = (shared_sub.T/shared_sub.sum(axis=1)).T.reset_index()

    palette = Dark2[8] + Set1[9] + Set2[8] + Set3[12] + Category20[20]

    offset = 200
    p = figure(plot_height=2*offset+max(500, shared_sub.shape[0]*50),
               plot_width=2*offset+1000,
               x_range=[0, 1.5], y_range=shared_sub.Group,
               min_border=offset,
               title="{} relative abundance across samples".format(rank),
               tooltips=list(zip(shared_sub.columns[1:], '@'+shared_sub.columns[1:]+'{0.00%}')))

    p.hbar_stack(shared_sub.columns[1:], y='Group', source=shared_sub,
                 height=.5, color=palette[:shared_sub.shape[1]-1],
                 legend_label=shared_sub.columns[1:].tolist())

    p.xaxis.visible = False
    p.grid.grid_line_color = None
    p.axis.major_label_text_font_size = "10pt"
    p.title.text_font_size = '14pt'

    output_file(outname)
    save(p)


def preproc_for_heatmap(shared, tax, meta=None, top=100):

    (data, _) = group_otus(shared, tax, 'Genus', top=top)

    # z-score first and biclustering
    z_data = (data - data.mean()) / data.std()

    # hierarchical clustering on both rows and columns
    try:
        row_links = leaves_list(linkage(z_data, method='average', metric='braycurtis', optimal_ordering=True))
    except ValueError:
        print('Something went wrong with the hierarchical clustering on rows')
        row_links = np.arange(z_data.shape[0])
    try:
        col_links = leaves_list(linkage(z_data.T, method='average', metric='braycurtis', optimal_ordering=True))
    except ValueError:
        print('Something went wrong with the hierarchical clustering on cols')
        col_links = np.arange(z_data.shape[1])

    abd_data = data.reset_index().melt(id_vars='Group').rename(columns={'value': 'abundance'})
    
    z_data = (z_data.iloc[row_links, col_links].T
              .reset_index()
              .melt(id_vars=['Genus'])
              .merge(abd_data)
              .rename(columns={'value': 'z_score', 'Group': 'sampleID', 'Genus': 'otu'})
              .set_index('sampleID')
    )

    info = {'otu': tax.groupby('Genus').agg('first')}

    if meta is not None:
        info['sample'] = meta[np.isin(meta.index, data.sampleID.unique())]

    return (z_data, info)

def clustermap(data, info, fig_dir='.', title=None, thresh=100, otu_groups='otu'):

    data = data.merge(info['otu'], left_on='otu', right_index=True)
    top = len(data.otu.unique())

    tooltips = dict(zip(info['otu'].columns, '@' + info['otu'].columns))
    tooltips[otu_groups] = '@otu'

    if 'sample' in info:
        data = data.merge(info['sample'], left_index=True, right_index=True)
        sample_tips = dict(zip(info['sample'].columns, '@' + info['sample'].columns))
        tooltips.update(sample_tips)

    cols = data.drop(['otu'], axis=1).columns
    other_tips = dict(zip(cols, '@' + cols))
    tooltips.update(other_tips)

    if title is None:
        title = "biclutered_heatmap_sample-by-genus-{}_top-{}".format(thresh, top)

    data.index.name = 'index'

    p = figure(plot_height=max(500, 10*len(data.index.unique())), plot_width=max(500, top*10),
               title=title,
               x_range=data.otu.unique(), y_range=data.index.unique().tolist(),
               x_axis_location="above",
               tools=TOOLS, tooltips=list(tooltips.items()))

    mapper = LinearColorMapper(palette='Viridis256',
                               low=data.z_score.min(),
                               high=data.z_score.max())

    p.rect(x="otu", y="index", width=1, height=1, source=data.reset_index(),
           line_color=None, fill_color=transform('z_score', mapper))

    color_bar = ColorBar(color_mapper=mapper, formatter=p.xaxis.formatter,
                         ticker=BasicTicker(desired_num_ticks=5),
                         location=(0,0))
    
    p.add_layout(color_bar, 'right')

    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "8pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = 'vertical'
    p.title.text_font_size = '14pt'

    output_file("{}/{}.html".format(fig_dir, title))
    save(p)
    
if __name__ == '__main__':
    args = parse_args()

    data = data(args)

    scatter(rank='Phylum', thresh=args.thresh, **data)
    scatter(rank='Class', thresh=args.thresh, **data)    
    stackplot(rank='Phylum', thresh=args.thresh, **data)
    stackplot(rank='Class', thresh=args.thresh, **data)

    (data, info) = preproc_for_heatmap(**data, top=args.ntop)
    clustermap(data, info, thresh=args.thresh, otu_groups='Genus')
    
