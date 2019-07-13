import pandas as pd
from glob import glob
import subprocess
from os.path import basename,dirname,splitext
import re

def guess_file_extension(files):
    extensions = { f.split(".")[-1] for f in files }

    if 'shared' in extensions:
        return ".shared"
    else:
        return ".count_table"
    
def get_file_extension(filename):
    splits = filename.split(".")

    if splits[-1]=='gz':
        return "{}.{}".format(*splits[-2:])
    return splits[-1]

class SequenceCounter:

    def __init__(self, name=None, path=None):
        self.name = name
        self.path = path
        self.extension = None
        self.fillInfo()
        self.setFileType()
        
    def fillInfo(self):
        if self.name is None:
            self.name = basename(glob(dirname(self.path))[0])
        self.step_nb = float(self.name.split('-')[0])
        try:
            self.id_threshold = re.findall(r"\d+", self.pattern)[0]
        except:
            pass

    def setFileType(self):
        files = glob(self.path)

        # Special case with taxa filter where the extension can be .count_table XOR shared
        if self.path.split(".")[-1] == "*":
            self.extension = guess_file_extension(files)
            self.compressed = False
            self.path += self.extension[1:]                
            return
        
        if len(files)==0:
            file_no_ext,ext = splitext(self.path)
        else:
            file_no_ext,ext = splitext(glob(self.path)[0])

        if ext == '.gz':
            self.compressed = True
            self.extension = splitext(file_no_ext)[1]
        else:
            self.compressed = False
            self.extension = ext

    def run(self,id_threshold,sample_names=None):
        files = glob(self.path)
        print(self.name, "{} files".format(len(files)))

        if self.extension.startswith(".fast"):
            counts = dict([ (basename(fastx).split("_R1")[0].replace("-","_"),
                             self.countFastx(fastx))
                            for fastx in files ])

            if sample_names is None:
                sample_names = list(counts.keys())
            summary_series = pd.Series(counts,name=self.name,index=sample_names).fillna(0).astype(int)

            return summary_series
        
        elif self.extension in [".shared",".count_table",".csv",".tsv"]:
            counts = [ self.countTable(table) for table in files ]
            if len(counts) > 1:
                counts = [ count for count in counts if id_threshold in count.name ]
            summary_df = pd.concat(counts,axis=1).reindex(sample_names).fillna(0)
                
            return summary_df 
        else:
            return pd.Series([])
            
    def countFastx(self,filename):
        open_cmd = "zcat" * (self.compressed) + "cat" * (not self.compressed)
        open_process = subprocess.Popen((open_cmd, filename), stdout=subprocess.PIPE)
        
        result = int(subprocess.check_output(('wc','-l'), stdin=open_process.stdout)
                     .decode()
                     .strip())
        if self.extension == ".fastq":
            return int(result / 4)
        elif self.extension == ".fasta":
            return int(result / 2)
        else:
            print("Wrong extension (neither .fasta or .fastq): {}".format(self.extension))

    def countTable(self,filename):
        name = self.name

        try:
            id_threshold = float(re.findall(r"\d[.]*[\d]*",
                                            splitext(basename(filename))[0])[0])
            name = "{0}_{1:d}".format(self.name,int(id_threshold))
        except:
            id_threshold = ""
                    
        if self.extension == ".shared":
            table = pd.read_csv(filename, index_col="Group", sep='\t').drop(["label","numOtus"],axis=1).T
        elif self.extension == ".count_table":
            table = pd.read_csv(filename, index_col=0, sep="\t").drop("total",axis=1)
        elif self.extension == ".csv":
            table = pd.read_csv(filename,index_col=0)
        elif self.extension == ".tsv":
            table = pd.read_csv(filename,index_col=0,sep='\t')
        else:
            print("Wrong extension (neither .shared or .count_table): {}".format(self.extension))

        summary = table.agg(lambda x: "{} ({} uniques)".format(x.sum(),(x>0).sum()))
        summary.index = summary.index.astype(str)
                
        return summary.rename(name)
