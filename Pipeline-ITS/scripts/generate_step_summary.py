from Bio.SeqIO.FastaIO import SimpleFastaParser

import pandas as pd
from glob import glob
import subprocess
import os
import re

def split_path(path):
    file_no_ext,ext = os.path.splitext(path)
    if ext == '.gz':
        file_no_ext, ext = os.path.splitext(file_no_ext)
        ext = "{}.gz".format(ext)
    
    path,fileRad = os.path.split(file_no_ext)

    return path,fileRad,ext

def get_thresh(filename):
    try:
        _,fileRad,_ = split_path(filename)
        thresh = re.findall(r"\d[.]*[\d]*", fileRad)[0]
    except:
        thresh = ""

    if thresh == "":
        thresh = "100"

    return thresh

def count_reads(filepath):
    path,fileRad,ext = split_path(filepath)
    name = re.split("_R\d",fileRad)[0] 
    filetype = "fast{}".format(ext.split(".")[1][-1]) # To handle the fasta/fa or fastq/fq

    cmd = "zgrep" * (ext.endswith(".gz")) + "grep" * (not ext.endswith(".gz"))
    pattern = "@M0" * (filetype=="fastq") + ">" * (filetype=="fasta")
    
    result = subprocess.run([cmd, '-c', pattern, filepath], stdout=subprocess.PIPE)
    count = int(result.stdout.decode().replace("\\n",""))
    return (name, count)

def count_unique_reads(filepath):
    path,fileRad,ext = split_path(filepath)
    name = re.split("_R\d",fileRad)[0]
    
    abundances = [ int(title.split("size=")[1].replace(";",""))
                   for (title,_) in SimpleFastaParser(open(filepath)) ]
    return (name, "{} ({} uniques)".format(sum(abundances),len(abundances)))
    

def process_mothur(step,filepath):
    _,fileRad,ext = split_path(filepath)

    if ext == ".shared":
        table = pd.read_table(filepath, index_col="Group").drop(["label","numOtus"],axis=1).T
    else:
        table = pd.read_table(filepath, index_col=0).drop("total",axis=1)

    name = "{}{}".format(step, get_thresh(filepath))
    summary = table.agg(lambda x: "{} ({} uniques)".format(x.sum(),(x>0).sum())).rename(name)
    
    return summary

def process_csv(filepath,step):

    name = "{}{}".format(step,get_thresh(filepath))    
    data = pd.read_csv(filepath,index_col=0)

    summary = data.agg(lambda x: "{} ({} uniques)".format(x.sum(),(x>0).sum())).rename(name)
    return summary
    

def count_samples(info,root_dir=".",samples=None):
    folder,pattern = info
    step_nb = int(folder.split("-")[0])
    
    if step_nb == 0:
        files = glob("{}/*R1*".format(os.path.dirname(pattern)))
    else:
        files = glob("{}/{}/{}".format(root_dir,folder,pattern))
    print(folder,"{} files".format(len(files)))

    _,_,ext = split_path(pattern)

    if ext.startswith(".fast") and step_nb < 5:
        seqs_per_sample = [ count_reads(filepath)
                            for filepath in files ]
        return pd.Series(dict(seqs_per_sample),name=folder)
    
    elif ext == ".fasta":
        seqs_per_sample = [ count_unique_reads(filepath)
                            for filepath in files ]
        return pd.Series(dict(seqs_per_sample),name=folder)
        
    elif ext == ".csv":
        csv_counts = [ process_csv(filepath,folder) for filepath in files ]
        return pd.concat(csv_counts,axis=1,sort=True)

def write_summary(root_dir,data_dir):
    
    steps = [("0-RawData", data_dir),
             ("1-preQC","*R1*.fastq.gz"),
             ("2-ITSxpress","*R1*.fastq"),
             ("3-NandSmallSeqsFiltering","*R1*.fastq"),
             ("4-qualityFiltering","*R1*.fastq.gz"),
             ("5-Dereplication","*.fasta"), # --> Generate summary with derep + abundance
             ("6-ChimeraRemoval","*.fasta"), # --> same function
             ("7-Denoising","*.fasta"),                          
             ("8-Clustering","*.csv"),
             ("9-LULU_correction","curated_table*.csv")                                         
    ]
    
    res_all_samples = [count_samples(step,root_dir=root_dir)
                       for step in steps]

    for res in res_all_samples:    
        res.index.name = "Sample"
        res.index = res.index.astype(str)    
    
    summary = pd.concat(res_all_samples,axis=1,sort=True)
    
    summary.index.name = "SampleID"

    summary.to_csv("sequences_per_sample_per_step.tsv",sep="\t")

if __name__ == '__main__':

    write_summary("/home/cedric/projects/HIDSI/1-Pipelines/2019-02-04-ITS-results",
                  "/home/cedric/projects/HIDSI/0-Data/ITS_testdata/*R1*.fastq.gz")
