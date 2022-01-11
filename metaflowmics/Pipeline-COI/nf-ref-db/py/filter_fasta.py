from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, default="iBOL_COI.fna")
    parser.add_argument('--tax', type=str, default="taxonomy_clean.tsv")
    args = parser.parse_args()

    return args

args = parse_args()

accessions = pd.read_csv(args.tax, sep="\t", index_col=0).index.astype(str)

with open(args.fasta, "r") as reader, open("rawSeq.fasta", "w") as writer:
    for (title, seq) in SimpleFastaParser(reader):
        acc = title.split()[0]

        if acc in accessions:
            writer.write(f">{acc}\n{seq}\n")
