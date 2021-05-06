import argparse

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
import h5py

NUCLS = list('ACGTN')
MAPPING_BASE5 = str.maketrans(''.join(NUCLS), '01234')
RC_TRANS = str.maketrans('ACGT', 'TGCA')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastqs', type=str, nargs='+')
    parser.add_argument('--split', type=int)
    parser.add_argument('--rc', action='store_true', default=False)
    args = parser.parse_args()

    return args

def seq2str(seq):
    if hasattr(seq, 'seq'):
        return str(seq.seq)
    return seq

def seq2intlist(bc):
    bc_vec = np.fromiter(seq2str(bc).translate(MAPPING_BASE5), dtype='uint8')
    return bc_vec

def seq2idx(bc):
    bc_idx = int(seq2str(bc).translate(MAPPING_BASE5), 5)
    return bc_idx

def get_lengths(filename):
    handle = open(filename)
    rid_len = len(next(handle).strip().split()[0])
    bc_len = len(next(handle).strip())

    return (rid_len, bc_len)

def fastqs_to_h5(fastqs, n_reads, bc_len=None, rid_len=None, split=0, rc=False):
    parsers = [FastqGeneralIterator(open(fastq)) for fastq in fastqs]

    data = [
        {'rid': np.empty(n_reads, dtype='S{}'.format(rid_len+10)),
         'index': np.zeros(n_reads, dtype='uint32'),
         'seq': np.zeros((n_reads, bc_len), dtype='uint8'),
         'qual': np.zeros((n_reads, bc_len), dtype='uint8')}
        for _ in fastqs]

    for k, parser in enumerate(parsers):
        for i, (rid, seq, qual) in enumerate(parser):
            if rc:
                seq = seq.translate(RC_TRANS)
            data[k]['rid'][i] = rid.split()[0]
            data[k]['index'][i] = seq2idx(seq)
            data[k]['seq'][i] = seq2intlist(seq)
            data[k]['qual'][i] = [ord(score)-33 for score in qual]

            if k == 1:
                (rid_fwd, rid_rev) = data[0]['rid'][i], data[1]['rid'][i]
                if rid_fwd != rid_rev:
                    print('Error: {} != {}: forward and reverse read IDs do not match.'
                          .format(rid_fwd, rid_rev))
                    exit(42)
        
    h5_handle = h5py.File('data_{}.h5'.format(split), 'w')

    options = ['fwd', 'rev'][:len(fastqs)]
    for k, orient in enumerate(options):
        for key, val in data[k].items():
            h5_handle.create_dataset('{}/{}'.format(key, orient), data=val)
    h5_handle.close()

def main():
    
    args = parse_args()
    n_reads = sum(1 for _ in open(args.fastqs[0])) // 4
    (rid_len, bc_len) = get_lengths(args.fastqs[0])
    
    fastqs_to_h5(args.fastqs, n_reads, bc_len=bc_len, rid_len=rid_len, split=args.split, rc=args.rc)

if __name__ == '__main__':
    main()
