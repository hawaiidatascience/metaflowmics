import pandas as pd
from os.path import basename,splitext
import re

class SequenceCounter:

    def __init__(self, name=None, path=None):
        self.name = name
        self.path = path
        self.extension = path.split('.')[-1]
        
    def run(self,sample_names=None):
        if self.extension in ["shared","count_table","csv","tsv"]:
            counts = self.countTable(self.path)
            summary_df = counts.reindex(sample_names).fillna(0)
            
            return summary_df 
        else:
            return pd.Series([])
            
    def countTable(self,filename):
        name = self.name

        try:
            id_threshold = float(re.findall(r"\d[.]*[\d]*",
                                            splitext(basename(filename))[0])[0])
            name = "{0}_{1:d}".format(self.name,int(id_threshold))
        except:
            id_threshold = ""
                    
        if self.extension == "shared":
            table = pd.read_csv(filename, sep='\t', dtype={'Group': str}).drop(["label","numOtus"],axis=1).set_index('Group').T
        elif self.extension == "count_table":
            table = pd.read_csv(filename, index_col=0, sep="\t").drop("total",axis=1)
        elif self.extension == "csv":
            table = pd.read_csv(filename,index_col=0)
        elif self.extension == "tsv":
            table = pd.read_csv(filename,index_col=0,sep='\t')
        else:
            print("Wrong extension (neither .shared or .count_table): {}".format(self.extension))
            
        summary = table.apply(lambda x: "{} ({} uniques)".format(x.sum(),(x>0).sum()))
        summary.index = summary.index.astype(str)
                
        return summary.rename(name)
