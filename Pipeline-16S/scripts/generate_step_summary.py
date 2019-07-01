import pandas as pd
from counter import SequenceCounter

def write_summary(root_dir,data_dir,id_threshold):
    
    steps = [("0-rawData",data_dir),
             ("*filterAndTrim","*R1*.fastq*"),
             ("*errorModel","*R1*.RDS"),
             ("*multipleSequenceAlignment","all_MSA.count_table"),
             ("*chimeraRemoval","all_chimera.count_table"),
             ("*preClassification","*.taxonomy"),
             ("*taxaFilter","all_taxaFilter*.*"),
             ("*subsampling","all_subsampling.count_table"),
             ("*clustering","all_clustering_*.shared"),
             ("*consensusClassification","*.taxonomy"),             
             ("*lulu","lulu_table_*.csv"),
             ("*singletonFilter","singleton_filtered_*.shared")
    ]
    
    denoising_step = pd.read_csv("count_summary.tsv", index_col=0, sep='\t')

    res_all_samples = sorted([
        SequenceCounter(name=name,path=pattern.replace("{1,2}","1")).run(id_threshold) if i==0 else
        SequenceCounter(path="{}/{}/{}".format(root_dir,name,pattern)).run(id_threshold)
        for i,(name,pattern) in enumerate(steps)
    ], key=lambda x: x[0])

    res_all_samples = list(map(lambda x: x[1], res_all_samples))

    res_all_samples.insert(3,denoising_step)

    for i,res in enumerate(res_all_samples):
        if res.size == 0:
            res_all_samples[i] = pd.Series(index=denoising_step.index,
                                           name=res.name)
            
        res.index.name = "Sample"

    summary = pd.concat(res_all_samples,axis=1,sort=True) #.fillna(method='bfill')
    summary.index.name = "SampleID"

    summary.iloc[:,-1] = [ '0' if pd.isnull(x) else x for x in summary.iloc[:,-1] ]
    summary = summary.fillna(method='bfill',axis=1)

    # for i in range(1,summary.shape[0]):
    #     for j in range(summary.shape[1]-1):
    #         if summary.iloc[i,j] == summary.iloc[i,j-1] and summary.iloc[i,j] != 0:
    #             summary.iloc[i,j] = '-'
    summary_clean = summary[ (summary.shift(1,axis=1) != summary) | (summary == '0') ].fillna('-')
    summary_clean.iloc[:,-1] = summary.iloc[:,-1]

    summary_clean.to_csv("sequences_per_sample_per_step_{}.tsv".format(id_threshold),sep="\t")
