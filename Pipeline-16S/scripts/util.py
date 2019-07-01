#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd

def getSubsamplingThreshold(countFile,quantile,minValue):
    table = pd.read_csv(countFile, index_col=0, sep="\t").drop("total",axis=1)
    sample_sizes = table.sum(axis=0).sort_values(ascending=False) 
    threshold = int(sample_sizes.quantile(quantile))

    if threshold > minValue:
        print(threshold)
    else:
        if sample_sizes[1] > minValue:
            print(minValue)
            ## print(minValue + " 1000") ==> for later
        elif sample_sizes[1] > 1000:
            print(1000)
        else:
            print(threshold)    

def filterAbundance(abundance,minAbundance=1):
    abundanceTable = pd.read_csv(abundance,index_col=0)
    abundanceTable = abundanceTable.loc[abundanceTable.sum(axis=1) > minAbundance]
    abundanceTable.to_csv("curated_{}".format(abundance))

def filterIds(abundance,fasta,taxonomy,threshold):
    ids = pd.read_csv(abundance,index_col=0).index.tolist()
    fasta = [ seq for seq in SeqIO.parse(fasta,"fasta") if seq.id in ids ]
    SeqIO.write(fasta,"singleton_filtered_{}.fasta".format(threshold),"fasta")

    taxonomyTable = pd.read_table(taxonomy,index_col=0)
    taxonomyTable.loc[ids].to_csv("singleton_filtered_{}.taxonomy".format(threshold),sep="\t")

def csvToShared(csvAbundance, threshold):
    shared = pd.read_csv(csvAbundance, index_col=0).T
    shared.index.name = "Group"
    shared.insert(loc=0, column="numOtus", value=[shared.shape[1]]*shared.shape[0])
    shared.reset_index(inplace=True)
    shared.insert(loc=0, column="label", value=[threshold]*shared.shape[0])
    shared.to_csv("singleton_filtered_{}.shared".format(threshold), sep="\t", index=False)


