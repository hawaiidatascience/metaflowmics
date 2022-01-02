import pandas as pd
from pathlib import tax
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tax', type=str)
args = parser.parse_args()

tax_file = Path(args.tax)

df = pd.read_csv(tax_file, sep="\t", header=None)
df.columns = ["acc", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
df = df.set_index('acc')

# remove duplicates (case doesn't matter)
duplicated = df.applymap(lambda x: x.lower()).duplicated()
df = df[~duplicated]

# def fix_tax(row):
#     unknowns = [i for i,v in enumerate(row.values) if "unknown" in v][::-1]
#     while unknowns:
#         idx = unknowns.pop()
#         prev = row.iloc[idx-1] # assume the highest level is always known
#         prev_rank = row.index[idx-1][0].lower()
#         prev_comps = prev.split("__")

#         new_name = f"{prev_rank}__{prev_comps[-1]}"
#         if len(prev_comps) > 1:
#             prefix = "__".join(prev_comps[:-1])
#             new_name = f"{prefix}__{new_name}"
#         row.loc[row.index[idx]] = new_name
#     return row

# df = df.apply(fix_tax, axis=1)
df.to_csv(f"{tax_file.parent}/{tax_file.stem}_clean.tsv", sep="\t")
