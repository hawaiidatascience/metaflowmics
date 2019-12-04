import argparse
from glob import glob

import pandas as pd
import numpy as np
import h5py
from sklearn.isotonic import IsotonicRegression

import matplotlib.pyplot as plt
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text')

NUCLS = list('ACGTN')
MAPPING_BASE5 = str.maketrans(''.join(NUCLS), '01234')
    
def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', type=str)
    parser.add_argument('--n_bases', type=int, default=int(1e6))
    parser.add_argument('--max_dist', type=int, default=2, help='Maximum hamming distance for alignments')
    args = parser.parse_args()

    return args

def seq2str(seq):
    if hasattr(seq, 'seq'):
        return str(seq.seq)
    return seq

def seq2intlist(bc):
    bc_vec = np.fromiter(seq2str(bc).translate(MAPPING_BASE5), dtype='uint8')
    return bc_vec

def get_h5_keys():
    h5_file = glob("*.h5")[0]
    handle = h5py.File(h5_file)

    return list(handle.get('seq').keys())

def load_h5():

    orientations = get_h5_keys()

    indexes = {orient: [] for orient in orientations}
    sequences = {orient: [] for orient in orientations}
    qualities = {orient: [] for orient in orientations}
    
    for h5 in sorted(glob("*.h5")):
        handle = h5py.File(h5)
        for orient in orientations:
            indexes[orient].append(handle.get('index/{}'.format(orient))[:])
            sequences[orient].append(handle.get('seq/{}'.format(orient))[:])            
            qualities[orient].append(handle.get('qual/{}'.format(orient))[:])
        handle.close()

    fastq_data = {
        orient: {'index': np.concatenate(indexes[orient]),
                 'seq': np.vstack(sequences[orient]),
                 'qual': np.vstack(qualities[orient])}
        for orient in orientations
    }
    return fastq_data

def idx2seq(nb):
    repr_b5 = "{:08}".format(int(np.base_repr(nb, 5)))
    nucl = repr_b5.translate(str.maketrans("01234", ''.join(NUCLS)))
    return nucl

def compute_transition_matrix(fastq_data, barcodes, max_dist=2, n_bases_max=1e5, eps=1e-6):
    '''
    - Barcodes in a numpy array of nucleotide index (base 5)
    - 
    '''

    # barcodes_index = np.array([int(''.join(map(str, bc)), 5) for bc in barcodes])
    bc_len = len(barcodes[0])
    transition_matrix = np.zeros((4, 5, 40)) # Transitions (A,C,G,T) to (A,C,G,T,N) and quality_scores
    generator = zip(fastq_data['index'], fastq_data['seq'], fastq_data['qual'])
    n_bases = 0

    mem = {}
    while n_bases < n_bases_max:

        if n_bases % 10*bc_len == 0: 
            print("{:,}/{:,}".format(n_bases, n_bases_max), end='\r')
        
        (code, intlist, qual) = next(generator)

        if code in mem:
            matching_bc = mem[code]
        else:
            matching_bc = barcodes[np.sum(barcodes != intlist[None, :], axis=1) <= max_dist]
            mem[code] = matching_bc

        if len(matching_bc) == 0:
            continue

        for bc in matching_bc:
            tup, counts = np.unique(np.vstack((bc, intlist, qual-1)),
                                    return_counts=True, axis=1)
            transition_matrix[tup[0, :], tup[1, :], tup[2, :]] += counts

            n_bases += bc_len

    sums = transition_matrix.sum(axis=0)
    sums[sums==0] = eps

    transition_matrix /= sums
    
    return transition_matrix

def regression_model(freqs, deg=2, same=False, method='isotonic'):
    '''
    - qual: all measured quality scores when a transition was observed
    - proportion of transitions for a given quality
    '''

    observed_transitions = (~np.isnan(freqs)) & (freqs>0)

    x = np.arange(1, 41)[observed_transitions]
    y = -np.log10(freqs[observed_transitions])

    if len(x) == 0:
        return np.zeros(40)

    if method == 'polynomial':
        z = np.polyfit(x, y, 3)
        polynom = np.poly1d(z)
        y_interp = 10**-polynom(np.arange(1, 41))
    elif method == 'isotonic':
        ir = IsotonicRegression(y_min=0, out_of_bounds='clip', increasing=not same)
        ir.fit(x, y)
        y_interp = 10**-ir.predict(np.arange(1, 41))
    else:
        print('Unknown method: {}. Aborting.'.format(method))
        exit(1)
    return y_interp

def compute_error_models(idx, barcodes, args, orient='fwd'):

    barcodes_int = np.array(
        [list(seq.translate(MAPPING_BASE5)) for seq in barcodes]
    ).astype('uint8')

    transitions = compute_transition_matrix(idx, barcodes_int, max_dist=args.max_dist, n_bases_max=args.n_bases)

    transitions_interp = np.array(
        [[regression_model(err_prob, same=(i==j))
          for i, err_prob in enumerate(transition_from_orig)]
         for j, transition_from_orig in enumerate(transitions)])

    display_models(transitions, transitions_interp, orient=orient)
    print('Error model computed')
    
    return transitions_interp

def display_models(trans_matrix, trans_matrix_interp, orient='fwd'):
    """
    """

    # Data stores
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

    fig, axes = plt.subplots(4, 5)
    for i, nucl_i in enumerate(NUCLS[:-1]):
        for j, nucl_j in enumerate(NUCLS):
            transition_label = "{}->{}".format(nucl_i, nucl_j)
            if transition_label not in data.index:
                continue
            
            ax = axes[i, j]
            data_s = data.loc[transition_label].sort_values(by='quality')
            
            ax.scatter(data_s.loc['y1', 'quality'],
                             data_s.loc['y1', 'P_err'], alpha=0.7, edgecolor='k')
            ax.plot(data_s.loc['y_hat', 'quality'],
                          data_s.loc['y_hat', 'P_err'], 'r-.')
            ax.set_title(r"{}->{}".format(nucl_i, nucl_j))

            if i < len(NUCLS)-2:
                ax.set_xticklabels([])

    plt.savefig("error_model_{}.pdf".format(orient), transparent=True)
    
def main():
    '''
    '''

    args = parse_args()

    barcodes = (pd.read_csv(args.meta, dtype=str, names=["sample_name", "fwd", "rev"])
                .dropna(axis=1, how='all')
                .drop('sample_name', axis=1))

    fastq_data = load_h5()

    transition_h5 = h5py.File('transition_probs.h5', 'w')
    for (orient, val) in fastq_data.items():
        transition_h5.create_dataset(orient, data=compute_error_models(val, barcodes[orient].unique(), args, orient=orient))
    transition_h5.close()
        
if __name__ == '__main__':
    main()
