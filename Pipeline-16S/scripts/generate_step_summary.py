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

def get_thresh(filename,removePath=True):
    try:
        if removePath:
            _,fileRad,_ = split_path(filename)
        else:
            fileRad = filename
        name = re.findall(r"\d[.]*[\d]*", fileRad)[0]
        return name
    except:
        return ""

def count_reads(filepath):
    path,fileRad,ext = split_path(filepath)
    name = re.split("_R\d",fileRad)[0] 
    filetype = "fast{}".format(ext.split(".")[1][-1]) # To handle the fasta/fa or fastq/fq

    cmd = "zgrep" * (ext.endswith(".gz")) + "grep" * (not ext.endswith(".gz"))
    pattern = "@M0" * (filetype=="fastq") + ">" * (filetype=="fasta")
    
    result = subprocess.run([cmd, '-c', pattern, filepath], stdout=subprocess.PIPE)
    count = int(result.stdout.decode().replace("\\n",""))
    return (name, count)

def process_mothur(step,filepath):
    _,fileRad,ext = split_path(filepath)

    if ext == ".shared":
        table = pd.read_table(filepath, index_col="Group").drop(["label","numOtus"],axis=1).T
    else:
        table = pd.read_table(filepath, index_col=0).drop("total",axis=1)

    name = "{}_{}".format(step, get_thresh(filepath))
    summary = table.agg(lambda x: "{} ({} uniques)".format(x.sum(),(x>0).sum())).rename(name)
    
    return summary

def process_lulu(filepath):

    name = "lulu{}".format(get_thresh(filepath))    
    data = pd.read_csv(filepath,index_col=0)

    summary = data.agg(lambda x: "{} ({} uniques)".format(x.sum(),(x>0).sum())).rename(name)
    return summary
    

def count_samples(info,root_dir=".",samples=None):
    folder,pattern = info

    if folder.startswith("0"):
        files = glob("{}/*R1*".format(os.path.dirname(pattern)))
    else:
        files = glob("{}/{}/{}".format(root_dir,folder,pattern))
    print(folder,"{} files".format(len(files)))

    ext = pattern.split('.')[-1].lower()
    
    if ext.startswith("fa"):
        seqs_per_sample = [ count_reads(filepath)
                            for filepath in files ]
        return pd.Series(dict(seqs_per_sample),name=folder)
        
    elif ext.startswith("count_table") or ext.startswith("shared"):
        res = [ process_mothur(folder,filepath) for filepath in files ]
        if len(res) > 1:
            res = sorted(res, key=lambda x: float(x.name.split("_")[-1]))
        return pd.concat(res, axis=1)

    elif 'lulu' in folder.lower():
        res = sorted([ process_lulu(filepath) for filepath in files ],
                     key=lambda x: float(get_thresh(x.name,removePath=False)))
        return pd.concat(res,axis=1)
    
    else:
        return pd.Series(["NA"]*len(samples), index=samples, name=folder)

def write_summary(root_dir,data_dir):
    
    steps = [("0-rawData",data_dir),
             ("1-filterAndTrim","*R1*.fastq*"),
             ("2-errorModel","*R1*.RDS"),
             ("5-multipleSequenceAlignment","all_MSA.count_table"),
             ("6-chimeraRemoval","all_chimera.count_table"),
             ("7-subsampling","all_subsampling.count_table"),
             ("8-clustering","all_clustering_*.shared"),
             ("9-consensusClassification","*.taxonomy"),             
             ("10-lulu","lulu_table_*.csv"),
             ("11-singletonFilter","abundance_*.shared"),                          
             ("12-taxaFilter","all_taxaFilter_*.shared"),             
    ]
    
    denoising_step = pd.read_table("count_summary.tsv", index_col=0)

    res_all_samples = [count_samples(step,root_dir=root_dir,samples=denoising_step.index)
                       for step in steps] 

    res_all_samples.insert(3,denoising_step)
    summary = pd.concat(res_all_samples,axis=1,sort=True)
    
    summary.index.name = "SampleID"

    summary.to_csv("sequences_per_sample_per_step.tsv",sep="\t")
