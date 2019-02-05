#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
from glob import glob

def removeNseqs(fastq,minLen,pairId):
    seqs = [ seq for seq in SeqIO.parse(fastq,"fastq")
             if "N" not in seq.seq or len(seq)<minLen ]
    SeqIO.write(seqs,pairId+"_noN.fastq","fastq")

def mergeSamplesFa(pattern="*.fasta"):
    '''
    Merges and rename fastas for all samples before making OTU table
    '''
    sequences = []
    for fasta_file in glob(pattern):
        sequences += [ seq for seq in SeqIO.parse(fasta_file,'fasta') ]
        
    for seq in sequences:
        seq.id = 'sample=' + seq.id
        seq.name = 'sample=' + seq.name
        
    SeqIO.write(sequences,'all_samples_merged.fasta','fasta')

    
def extractFastaLulu(fasta,idsToExtract,threshold):
    ids = pd.read_csv(idsToExtract,
                      header=None)[0].values
    sequences = [ seq for seq in SeqIO.parse(fasta,"fasta")
                  if seq.id.split(";")[0] in ids ]
    SeqIO.write(sequences,"lulu_fasta{}.fasta".format(threshold),"fasta")
    
def fillAbundances(abundanceTable,taxonomy,threshold):
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
