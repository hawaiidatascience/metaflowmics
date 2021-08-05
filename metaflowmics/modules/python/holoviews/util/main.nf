
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process HOLOVIEWS_PREPARE {
    tag "$otu_id"
    label "process_medium"
    label "plot"

    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::datatable pandas>=1" : null)

    input:
    tuple val(otu_id), file(shared), file(tax)

    output:
    tuple val(otu_id), file("metagenome.csv")

    script:
    def software = getSoftwareName(task.process)
    def max_mem = task.memory ? task.memory.getBytes() : "None"
    """
    #!/usr/bin/env python

    import datatable as dt
    import pandas as pd

    # Load taxonomy
    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    otu_meta = dt.fread("$tax", columns=dict(Size=None)).to_pandas()
    otu_meta = otu_meta.set_index("OTU").Taxonomy.str.rstrip(";").str.split(";", expand=True)
    otu_meta.columns = ranks[:otu_meta.shape[1]]

    # Load count data
    abund = dt.fread("$shared", columns=dict(label=None, numOtus=None), 
                     memory_limit=$max_mem)
    abund = abund.to_pandas().melt(id_vars="Group", var_name="OTU", value_name="count")

    # binary columns are set to bool by fread --> convert to int
    abund = abund.astype(dict(count=int))[abund["count"]>0]

    # compute relative abundance
    abund["relabund"] = abund.groupby("Group")["count"].transform(lambda x: x/sum(x) if any(x>0) else 0)
    
    # Merge
    data = abund.merge(otu_meta.reset_index(), on="OTU", how="left")
    
    data.to_csv("metagenome.csv", index=False)
    """
}
