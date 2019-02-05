#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
from glob import glob
import re

def mergeSamplesFa(pattern="*.fasta"):
    '''
    Merges and rename fastas for all samples before making OTU table
    '''
    sequences = []
    for fasta_file in glob(pattern):
        sequences += [ seq for seq in SeqIO.parse(fasta_file,'fasta') ]
        
    for seq in sequences:
        sample = re.split("_R\d",seq.id)[0]
        abundance = re.split("size=",seq.id)[-1]
        new_id = "sample={};size={}".format(sample,abundance)
        
        seq.id = new_id
        seq.name = new_id
        seq.title = ""
        
    SeqIO.write(sequences,'all_samples_merged.fasta','fasta')
    
def fillAbundances(abundanceTable,taxonomy,threshold):
    '''
    Fills taxonomy annotation in the abundance table
    '''
    annotations = pd.read_table(taxonomy,
                                header=None,
                                index_col=0)
    annotations.index = [ idx.split(';')[0]
                          for idx in annotations.index]

    table = pd.read_csv(abundanceTable,
                        index_col=0)
    taxonomy = annotations.loc[table.index]
    table["taxonomy"] = annotations.loc[table.index][1]

    table.to_csv("abundance_table_annotated_ID={}.tsv".format(threshold),
                 sep="\t")
