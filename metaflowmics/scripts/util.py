#!/usr/bin/env python
import pandas as pd

def get_subsampling_threshold(sharedFile, quantile, minValue, customValue):
    if customValue > 0:
        print(customValue)
        return

    table = (pd.read_csv(sharedFile, sep="\t", index_col=0)
             .drop(["total"],axis=1))
    
    sample_sizes = table.sum(axis=0).sort_values(ascending=False) 
    threshold = int(sample_sizes.quantile(quantile))

    if threshold > minValue:
        print(threshold)
    else:
        if sample_sizes[1] > minValue:
            print(minValue)
            ## TODO: enable multiple space delimited thresholds
        elif sample_sizes[1] > 1000:
            print(1000)
        else:
            print(threshold)    

