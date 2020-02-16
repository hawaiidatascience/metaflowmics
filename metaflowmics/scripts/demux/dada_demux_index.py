import argparse
from time import time

import pandas as pd
import numpy as np
import h5py
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

NUCLS = list('ACGTN')
MAPPING_BASE5 = str.maketrans(''.join(NUCLS), '01234')

def timer(func):
    def timed(*args, **kw):
        ts = time()
        result = func(*args, **kw)
        te = time()
        print('{0}  {1:2.2f}s'.format(func.__name__, (te - ts)))
        return result
    return timed
    
def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str)
    parser.add_argument('--error-model', type=str)
    parser.add_argument('--meta', type=str)
    parser.add_argument('--split', type=int, default=0)
    parser.add_argument('--max-mismatches', type=int, default=2)    
    args = parser.parse_args()

    return args

def seq2str(seq):
    if hasattr(seq, 'seq'):
        return str(seq.seq)
    return seq

def seq2intlist(bc):
    bc_vec = np.fromiter(seq2str(bc).translate(MAPPING_BASE5), dtype='uint8')
    return bc_vec

def idx2seq(nb):
    repr_b5 = "{:08}".format(int(np.base_repr(nb, 5)))
    nucl = repr_b5.translate(str.maketrans("01234", ''.join(NUCLS)))
    return nucl

def from_h5(filename, field=None):
    handle = h5py.File(filename, 'r')
    keys = handle.keys()

    if field is not None:
        keys = ["{}/{}".format(field, key) 
                for key in handle.get(field).keys()]
    
    data = {key.split('/')[-1]: np.array(handle.get(key)[:]) for key in keys}

    handle.close()

    return data

def calc_probs(sequences, qualities, error_model, barcodes, max_mismatches=2):
    '''
    data: index read data (seq index, codes and quals)
    '''

    print('Loading data')
    transition_probs = from_h5(error_model)
    orientations = transition_probs.keys()
    
    print('Computing likelihood')

    bc_sequences = {orient: np.array([seq2intlist(bc) for bc in barcodes[orient]])
                    for orient in orientations}

    transitions_all = [transition_probs[orient][
        bc_sequences[orient][:, None],
        sequences[orient],
        qualities[orient]-1].prod(axis=2) for orient in orientations]

    if len(transitions_all) == 1:
        assignments = transitions_all[0].argmax(axis=0)
    else:
        assignments = np.multiply(*transitions_all).argmax(axis=0)

    diffs = sum([bc_sequences[orient][assignments, :] != sequences[orient]
                 for orient in orientations]).sum(axis=1)

    assignments = barcodes.iloc[assignments].reset_index()
    assignments.insert(1, 'mismatches', diffs)
    assignments.loc[diffs > max_mismatches, 'sample_name'] = np.nan
    
    return assignments.values

def format_result(assignments, rids, seq):

    trans = str.maketrans('01234', 'ACGTN')
    index_string = np.array(
        [[''.join(x) for x in np.core.defchararray.translate(codes.astype('U1'), trans)]
        for codes in seq.values()])

    result = np.vstack([rids[None, :], index_string, assignments[:, [0, 1]].T]).T

    return pd.DataFrame(result)

def main():
    '''
    '''

    args = parse_args()

    barcodes = pd.read_csv(args.meta, dtype=str, names=["sample_name", "fwd", "rev"]).dropna(how='all', axis=1)
    barcodes.sample_name = barcodes.sample_name.str.replace("[^a-zA-Z0-9_.]", "_")
    barcodes.set_index('sample_name', inplace=True)

    read_ids = from_h5(args.data, field='rid')['fwd'].astype(str)
    sequences = from_h5(args.data, field='seq')
    quals = from_h5(args.data, field='qual')

    assignments = calc_probs(sequences, quals, args.error_model, barcodes, max_mismatches=args.max_mismatches)
    result = format_result(assignments, read_ids, sequences)

    result.to_csv("demux_info_{}.tsv".format(args.split), sep='\t', header=False, index=False)

    sample_sizes = result[result.columns[-2]].astype(str).value_counts()

    (
        sample_sizes[~sample_sizes.index.isnull()]
        .to_csv("sample_counts_{}.csv".format(args.split), header=False)
    )

if __name__ == '__main__':
    main()
