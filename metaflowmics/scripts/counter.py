from pathlib import Path
import re

import pandas as pd

class SequenceCounter:

    def __init__(self, name=None, path=None):
        self.name = name
        self.path = Path(path)
        
        otu_threshold = re.findall(r"[1-9][0-9]{1,2}", self.path.stem)

        if otu_threshold:
            self.name = "{}_{}".format(name, otu_threshold[0])
        
    def run(self, sample_names=None):
        if self.path.suffix in [".shared", ".count_table", ".csv", ".tsv"]:
            counts = self.countTable(self.path)
            summary_df = counts.reindex(sample_names).fillna(0)
            
        elif self.path.suffix == ".summary":
            summary_df = self.getSummary(self.path)
        else:
            summary_df = pd.Series([])

        return summary_df.rename(self.name)
            
    def countTable(self, filename):                    
        if self.path.suffix == ".shared":
            header = open(filename).readline().split('t')
            dtype = dict([(x, str) if x == 'Group' else (x, int) for x in header])

            table = (pd.read_csv(filename, sep='\t', dtype=dtype,
                                 keep_default_na=False, low_memory=False)
                     .drop(["label","numOtus"],axis=1)
                     .set_index('Group').T)
        elif self.path.suffix == ".count_table":
            table = pd.read_csv(filename, index_col=0, sep="\t").drop("total", axis=1)
        elif self.path.suffix == ".csv":
            table = pd.read_csv(filename, index_col=0)
        elif self.path.suffix == ".tsv":
            table = pd.read_csv(filename, index_col=0, sep='\t')
        else:
            print("Wrong extension (neither .shared or .count_table): {}".format(self.path.suffix))
            exit(42)
            
        summary = table.apply(lambda x: "{} ({} uniques)".format(x.sum(), (x > 0).sum()))

        return summary

    def getSummary(self, path, idx='group', cols=['nseqs', 'sobs']):
        
        summary = (pd.read_csv(path, sep='\t', usecols=[idx]+cols, dtype={idx: str})
                   .set_index(idx)
                   .astype(int)
                   .apply(lambda x: "{} ({} uniques)".format(x.nseqs, x.sobs), axis=1))
        return summary
