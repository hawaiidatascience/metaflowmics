#!/usr/bin/env python

import sys
from glob import glob
import pandas as pd

thresh = sys.argv[2]
shared = pd.read_table(sys.argv[1],index_col=0)
otus = shared.columns[2:]

accnos = glob("*.accnos")[0]
good_indices = [ f.strip() for f in open(accnos) ][1:]
cols = shared.columns[:2].tolist() + good_indices

shared.loc[:,'numOtus'] = [ len(good_indices) ] * shared.shape[0]
shared[cols].to_csv("all_taxaFilter_{}.shared".format(thresh),sep="\t")
