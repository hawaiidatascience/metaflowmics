from Bio.SeqIO.FastaIO import SimpleFastaParser

import argparse


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--query', type=str)
    parser.add_argument('--ref', type=str)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()

    return args

def sync(query, ref, output):

    ref_ids = [title.split()[0] for (title, _) in
               SimpleFastaParser(open(ref))]

    query_seq = {title.split()[0]: f">{title}\n{seq}\n" for (title, seq) in
                 SimpleFastaParser(open(query))}

    with open(output, "w") as writer:
        for prot_id in ref_ids:
            writer.write(query_seq[prot_id])        

if __name__ == '__main__':
    args = parse_args()

    sync(args.query, args.ref, args.output)

