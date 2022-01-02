import pandas as pd

def get_taxa(df, level=1, max_ct=1000, taxa=[]):

    if df.empty:
        return taxa

    lineage = df[level]
    if level > 1:
        lineage = df[level-1] + ';' + lineage

    freqs = lineage.value_counts()
    taxa.append(freqs[freqs<max_ct])
    large_groups = freqs.index[freqs>=max_ct]
    selection = lineage.isin(large_groups).to_numpy()
    df_lg = df[selection]

    print(f"level={level}, max_freq={freqs.max()}")

    print(freqs[freqs>=max_ct])

    import ipdb;ipdb.set_trace()
    return get_taxa(df_lg, level+1, max_ct, taxa)


info = pd.read_csv("iBOL-database/taxa.csv", header=None)
taxa = get_taxa(info, level=1, max_ct=50000)

print("all_taxa")
print(taxa)

all_taxa = [x.split(';')[-1][3:] for tax in taxa for x in tax.index]

with open("iBOL-database/taxa_to_dl.txt", "w") as writer:
    for t in all_taxa:
        if "undef" not in t:
            writer.write(f"{t}\n")
    
