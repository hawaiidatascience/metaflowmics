import argparse

import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.palettes import Set2, Set3, Category20
from bokeh.models import ColorBar, LinearColorMapper, BasicTicker
from bokeh.transform import transform

RANKS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
TOOLS = ['hover', 'box_zoom']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--shared', type=str)
    parser.add_argument('-t', '--taxonomy', type=str)
    parser.add_argument('--thresh', type=str)
    
    args = parser.parse_args()

    return args

def load_shared(path):
    df = pd.read_csv(path, index_col='Group', sep='\t').drop(['label', 'numOtus'], axis=1)
    return df

def load_tax(path):
    df = (pd.read_csv(path, index_col='OTU', sep='\t').drop('Size', axis=1).Taxonomy
          .str.strip(';')
          .str.replace('-', '_', regex=False)
          .str.split(';', expand=True))
    df.columns = RANKS[:df.shape[1]]

    return df

def data(args):

    return {'tax': load_tax(args.taxonomy),
            'shared': load_shared(args.shared)}

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
        p = figure(title=name.capitalize(), tooltips=tooltips)
        p.circle(x="prevalence", y="abundance", size=5, alpha=0.5, hover_color="red", source=data.loc[otus])
        plots.append(p)

    grid = gridplot(plots, ncols=5, plot_width=250, plot_height=250)
    save(grid)


def stackplot(shared, tax, thresh=100, rank='Phylum', top=20):

    if rank == 'Class':
        shared_sub = shared.loc[:, tax.Phylum == 'Proteobacteria']
        tax_sub = tax.loc[shared_sub.columns]
    else:
        shared_sub = shared.copy()
        tax_sub = tax.copy()

    shared_sub = shared_sub.T.groupby(tax_sub[rank]).agg('sum').T
    top = min(top, shared_sub.shape[1])

    prevalence = (shared_sub > 0).sum()
    abundance = shared_sub.sum()
    criteria = (prevalence * abundance).sort_values(ascending=False)

    shared_sub = shared_sub.loc[:, criteria.index[:top]]
    shared_sub = (shared_sub.T/shared_sub.sum(axis=1)).T.reset_index()

    palette = Set2[8] + Set3[12] + Category20[20]

    p = figure(plot_height=max(500, shared_sub.shape[0]*50), plot_width=1000,
               x_range=[0, 1.5], y_range=shared_sub.Group,
               title="{} relative abundance across samples".format(rank),
               tooltips=list(zip(shared_sub.columns[1:], '@'+shared_sub.columns[1:]+'{0.00%}')))

    p.hbar_stack(shared_sub.columns[1:], y='Group', source=shared_sub,
                 height=.5, color=palette[:top],
                 legend_label=shared_sub.columns[1:].tolist())

    output_file("barplot_OTU-{}_by-{}_top-{}.html".format(thresh, rank, top))
    save(p)

def clustermap(shared, tax, thresh=100, top=100):

    tax_info = tax.groupby('Genus').agg('first')
    prevalence = (shared > 0).sum()
    abundance = shared.sum()
    criteria = (prevalence * abundance).sort_values(ascending=False)

    data = shared.loc[:, criteria.index[:top]]
    data = data.T.groupby(tax.Genus).agg(sum)

    # z-score first and biclustering
    data = ((data.T - data.mean(axis=1)) / data.std(axis=1)).T

    # hierarchical clustering on both rows and columns
    try:
        row_links = linkage(data, method='average', metric='braycurtis', optimal_ordering=True)
        col_links = linkage(data.T, method='average', metric='braycurtis', optimal_ordering=True)
    except:
        return

    data = (
        data.iloc[leaves_list(row_links), leaves_list(col_links)]
        .reset_index()
        .melt(id_vars='Genus')
        .rename(columns={'value': 'z_score', 'Group': 'Sample'})
    )

    title = "biclutered_heatmap_OTU{}_top-{}".format(thresh, top)

    data = pd.concat([data.set_index('Genus'), tax_info.loc[data.Genus]], axis=1).reset_index()
    tooltips = [('Sample', '@Sample')] + list(zip(tax.columns, '@'+tax.columns)) + [('z-score', '@z_score')]

    p = figure(plot_width=max(500, shared.shape[0]*10), plot_height=max(500, top*10), title=title,
               x_range=data.Sample.unique(), y_range=data.Genus.unique(),
               x_axis_location="above",
               tools=TOOLS, tooltips=tooltips)

    mapper = LinearColorMapper(palette='Viridis256', low=data.z_score.min(), high=data.z_score.max())

    p.rect(x="Sample", y="Genus", width=1, height=1, source=data,
           line_color=None, fill_color=transform('z_score', mapper))

    color_bar = ColorBar(color_mapper=mapper, formatter=p.xaxis.formatter,
                         ticker=BasicTicker(desired_num_ticks=5),
                         location=(0,0))
    
    p.add_layout(color_bar, 'right')

    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "6pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = 'vertical'

    output_file("{}.html".format(title))
    save(p)
    
if __name__ == '__main__':
    args = parse_args()

    data = data(args)

    scatter(rank='Phylum', thresh=args.thresh, **data)
    scatter(rank='Class', thresh=args.thresh, **data)    
    stackplot(rank='Phylum', thresh=args.thresh, **data)
    stackplot(rank='Class', thresh=args.thresh, **data)    
    clustermap(**data, thresh=args.thresh)
    
