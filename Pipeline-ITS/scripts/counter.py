import pandas as pd
from glob import glob
import subprocess
from os.path import basename,dirname,splitext
import re

class SequenceCounter:

    def __init__(self,name,
                 folder=None,
                 pattern=None,
                 path=None):
        self.name = name
        self.folder = folder
        self.pattern = pattern
        self.path = path
        self.fillInfo()
        self.setFileType()

    def fillInfo(self):
        if self.path is None:
            self.path = "{}/{}".format(self.folder,self.pattern)
        elif self.pattern is None:
            self.folder = dirname(self.path)
            self.pattern = basename(self.path)
        else:
            print("Error: You need to provide a path or [folder+pattern]")
            exit(1)
        try:
            self.id_threshold = re.findall(r"\d[.]*[\d]*", self.pattern)[0]
        except:
            pass

    def setFileType(self):
        file_no_ext,ext = splitext(glob(self.path)[0])

        if ext == '.gz':
            self.compressed = True
            self.extension = splitext(file_no_ext)[1]
        else:
            self.compressed = False
            self.extension = ext

    def run(self):
        files = glob(self.path)
        print(self.name, "{} files".format(len(files)))
        
        if self.extension.startswith(".fast"):
            counts = [ (basename(fastx).split("_R1")[0],
                        self.countFastx(fastx))
                       for fastx in files ]
            return pd.Series(dict(counts),name=self.name)
        
        elif self.extension in [".shared",".count_table",".csv",".tsv"]:
            counts = [ self.countTable(table) for table in files ]
            if len(counts) > 1:
                counts = sorted(counts, key=lambda x: float(x.name.split("_")[-1]))
            return pd.concat(counts,axis=1)
            
    def countFastx(self,filename):
        open_cmd = "zcat" * (self.compressed) + "cat" * (not self.compressed)
        open_process = subprocess.Popen((open_cmd, filename), stdout=subprocess.PIPE)
        
        result = int(subprocess.check_output(('wc','-l'), stdin=open_process.stdout)
                     .decode()
                     .strip())
        if self.extension == ".fastq":
            return result / 4
        elif self.extension == ".fasta":
            return result / 2
        else:
            print("Wrong extension (neither .fasta or .fastq): {}".format(self.extension))

    def countTable(self,filename):
        name = self.name

        try:
            id_threshold = float(re.findall(r"\d[.]*[\d]*",
                                            splitext(basename(filename))[0])[0])
            name = "{0}_{1:.2g}".format(self.name,id_threshold)
        except:
            id_threshold = ""
                    
        if self.extension == ".shared":
            table = pd.read_table(filename, index_col="Group").drop(["label","numOtus"],axis=1).T
        elif self.extension == ".count_table":
            table = pd.read_table(filename, index_col=0).drop("total",axis=1)
        elif self.extension == ".csv":
            table = pd.read_csv(filename,index_col=0)
        elif self.extension == ".tsv":
            table = pd.read_table(filename,index_col=0)
        else:
            print("Wrong extension (neither .shared or .count_table): {}".format(self.extension))
                
        summary = table.agg(lambda x: "{} ({} uniques)".format(x.sum(),(x>0).sum()))

        return summary.rename(name)
