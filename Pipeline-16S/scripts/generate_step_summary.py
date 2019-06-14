import pandas as pd
from counter import SequenceCounter

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
             ("12-taxaFilter","abundance_table_*.shared"),             
    ]
    
    denoising_step = pd.read_table("count_summary.tsv", index_col=0)

    res_all_samples = [SequenceCounter(name,folder="{}/{}".format(root_dir,name),
                                       pattern=pattern).run() if i>0 else
                       SequenceCounter(name,path=data_dir.replace("{1,2}","1")).run()
                       for i,(name,pattern) in enumerate(steps)]

    res_all_samples.insert(3,denoising_step)

    for i,res in enumerate(res_all_samples):
        if res is None:
            res = pd.Series(index=denoising_step.index,
                            name=steps[i][0])
        elif res.size ==0:
            res_all_samples[i] = pd.Series(index=denoising_step.index,
                                           name=res.name)
            
        res.index.name = "Sample"
        res.index = [ row.split("_")[0] for row in res.index.astype(str) ]

    summary = pd.concat(res_all_samples,axis=1,sort=True)    
    summary.index.name = "SampleID"

    summary.to_csv("sequences_per_sample_per_step.tsv",sep="\t")
