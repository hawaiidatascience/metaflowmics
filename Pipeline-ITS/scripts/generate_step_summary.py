import pandas as pd
from counter import SequenceCounter

def write_summary(root_dir,data_dir):
    
    steps = [("0-RawData", data_dir),
             ("1-preFiltering","*R1*.fastq*"),
             ("2-ITSxpress","*R1*.fastq"),
             ("3-NandSmallSeqsFiltering","*R1*.fastq"),
             ("4-qualityFiltering","*R1*.fastq.gz"),
             ("8-Clustering","*.csv"),
             ("9-LULU_correction","curated_table*.csv")                                         
    ]

    denoising_step = pd.read_table("count_summary.tsv", index_col=0)
    
    res_all_samples = [SequenceCounter(name,folder="{}/{}".format(root_dir,name),
                                       pattern=pattern).run() if i>0 else
                       SequenceCounter(name,path=data_dir.replace("{1,2}","1")).run()
                       for i,(name,pattern) in enumerate(steps)]

    res_all_samples.insert(5,denoising_step)    
    
    for i,res in enumerate(res_all_samples):
        if res is None:
            res = pd.Series(index=denoising_step.index,
                            name=steps[i][0])
        res.index.name = "Sample"
        res.index = res.index.astype(str)

    summary = pd.concat(res_all_samples,axis=1,sort=True)
    
    summary.index.name = "SampleID"

    summary.to_csv("sequences_per_sample_per_step.tsv",sep="\t")
