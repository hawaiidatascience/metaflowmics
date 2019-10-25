import sys
import pandas as pd
from pathlib import Path

def removeEmptySamples(count_table_path):
    count_table = (pd.read_csv(count_table_path,
                               sep='\t',
                               dtype={'Representative_Sequence': str})
                   .set_index('Representative_Sequence'))
    nz_samples = count_table.drop('total', axis=1).values.sum(axis=0) > 0

    count_table = count_table.loc[:, [True] + nz_samples.tolist()]

    Path(count_table_path).unlink()
    count_table.to_csv(count_table_path, sep='\t')
    
if __name__ == '__main__':
    removeEmptySamples(sys.argv[1])
