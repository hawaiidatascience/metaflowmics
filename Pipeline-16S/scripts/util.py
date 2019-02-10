#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
    
def filterIds(ids,fasta,taxonomy,threshold):
    ids = [ name.strip() for name in open(ids,"r").readlines() ]
    fasta = [ seq for seq in SeqIO.parse(fasta,"fasta") if seq.id in ids ]
    SeqIO.write(fasta,"OTU_{}.fasta".format(threshold),"fasta")

    taxonomyTable = pd.read_table(taxonomy,index_col=0)
    taxonomyTable.loc[ids].to_csv("consensus_{}.taxonomy".format(threshold),sep="\t")

def filterAbundance(abundance,minAbundance=1):
    abundanceTable = pd.read_csv(abundance,index_col=0)
    abundanceTable = abundanceTable.loc[abundanceTable.sum(axis=1) > minAbundance]
    abundanceTable.to_csv("curated_{}".format(abundance))

def csvToShared(csvAbundance, threshold):
    shared = pd.read_csv(csvAbundance, index_col=0).T
    shared.index.name = "Group"
    shared.insert(loc=0, column="numOtus", value=[shared.shape[1]]*shared.shape[0])
    shared.reset_index(inplace=True)
    shared.insert(loc=0, column="label", value=[threshold]*shared.shape[0])
    shared.to_csv("abundance_{}.shared".format(threshold), sep="\t", index=False)
