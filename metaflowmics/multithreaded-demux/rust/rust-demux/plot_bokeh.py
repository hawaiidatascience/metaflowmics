import numpy as np
import pandas as pd

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot


NUCLS = list('ACGTN')


def display_models_bokeh(trans_matrix, trans_matrix_interp):
    data = []
    for matrix in [trans_matrix, trans_matrix_interp]:
        # For each matrix, stores non-zero probabilities with transitions and quality
        nz_info = np.nonzero(matrix)
        data.append(np.vstack((matrix[nz_info], *nz_info)).T)
    labels = np.array(['y1']*len(data[0]) + ['y_hat']*len(data[1]))
    data = np.vstack(data)
    transitions = ["{}->{}".format(NUCLS[i], NUCLS[j]) for i, j in data[:, [1, 2]].astype(int)]
    data = pd.DataFrame({'transition': transitions,
                         'quality': 1 + data[:, 3].astype(int),
                         'P_err': data[:, 0],
                         'type': labels})
    data.set_index(['transition', 'type'], inplace=True)
    plots = []
    tools = ['hover', 'box_zoom', 'reset']
    observed_transitions = set(pd.MultiIndex.get_level_values(data.index,'transition'))
    for i, nucl_i in enumerate(NUCLS):
        for j, nucl_j in enumerate(NUCLS[::-1]):
            transition_label = "{}->{}".format(nucl_i, nucl_j)
            if transition_label not in observed_transitions:
                plots.append(None)
                continue
            data_s = data.loc[transition_label].sort_values(by='quality')

            if "y_hat" not in data_s.index:
                plots.append(None)
                continue
                
            p = figure(title="{}->{}".format(nucl_i, nucl_j), tooltips=[], tools=tools,
                       x_range=(0, 40), y_range=(-0.1, 1.1))
            p.circle(x='quality', y='P_err', source=data_s.loc[['y1']], alpha=0.7, size=8)
            p.line(x='quality', y='P_err', source=data_s.loc[['y_hat']], color='red',
                   line_dash="dashed", line_width=2)
            plots.append(p)
    grid = gridplot(plots, ncols=5, plot_width=300, plot_height=300)
    output_file("error_model.html")
    save(grid)

    return data
        
