import argparse
from pathlib import Path
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


FIXES = dict(
    immidjanzen01="immidJanzen01",
    Enchytraeidaegen="EnchytraeidaeGEN",
    oecophor01="Oecophor01"
)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fna', type=str)
    parser.add_argument('--outdir', type=str, default="rdp")
    args = parser.parse_args()
    args.fna = Path(args.fna)

    return args

def make_taxonomy(fasta):
    tax = dict()

    lineages = set()
    with open(fasta) as reader:
        for (title, _) in SimpleFastaParser(reader):
            for (kw, repl) in FIXES.items():
                if kw in title:
                    title = title.replace(kw, repl)
            (acc, lin) = title.split()

            if lin.lower() not in lineages:
                tax[acc] = lin.split(";")
                lineages.add(lin.lower())

    tax_df = pd.DataFrame(tax).T
    tax_df.index.name = "acc"
    tax_df.columns = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    return tax_df

def filter_fasta(fasta, accessions, output):
    with open(fasta) as reader, open(output, "w") as writer:
        for (title, seq) in SimpleFastaParser(reader):
            acc = title.split()[0]
            if acc in accessions:
                writer.write(f">{acc}\n{seq}\n")

args = parse_args()

print('Preparing taxonomy')
tax = make_taxonomy(args.fna)
tax.to_csv(f"{args.outdir}/taxonomy.tsv", sep="\t")

print('Filtering sequences')
filter_fasta(args.fna, tax.index, output=f"{args.outdir}/rawSeq.fasta")
