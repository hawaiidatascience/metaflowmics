
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)


process HOLOVIEWS_CLUSTERMAP {
    tag "$otu_id"
    label "process_low"
    label "plot"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        meta:meta) }

    container "nakor/metaflowmics-python:0.0.1"
    conda (params.enable_conda ? "conda-forge::datatable pandas>=1" : null)

    input:
    tuple val(otu_id), file(mg)

    output:
    tuple val(otu_id), file("clustermap*.html")

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env python

    import pandas as pd
    from scipy.cluster.hierarchy import linkage, leaves_list
    import holoviews as hv

    hv.extension("bokeh")

    data = pd.read_csv("${mg}")
    
    # Filter OTUs
    scores = data.groupby("OTU").relabund.agg(max_abd=max, sd="std")
    otus = scores.index[(scores.max_abd > $params.min_abund) & (scores.sd > 0)]
    data = data[data.OTU.isin(otus)].set_index(["Group", "OTU"])

    # Compute z_scores
    zscores = data.relabund.unstack(fill_value=0)
    zscores = (zscores - zscores.mean()) / zscores.std()

    # Cluster rows and columns
    group_order = sorted(zscores.index)
    feature_order = sorted(zscores.columns)

    try:
        sample_links = leaves_list(linkage(
            zscores, method="average", metric="braycurtis"
        ))
        group_order = zscores.index[sample_links]
    except ValueError:
        print("Something went wrong with the sample clustering")

    try:
        feature_links = leaves_list(linkage(
            zscores.T, method="average", metric="braycurtis"
        ))
        feature_order = zscores.columns[feature_links]
    except ValueError:
        print("Something went wrong with the OTU clustering")
    
    # Reformat data
    data = data.merge(
        zscores.stack().rename("z_score"), 
        left_index=True, right_index=True, how="outer"
    ).reset_index()

    # Plot
    hm_opt = dict(
        tools=["hover"],
        height=max(10*data.Group.nunique(), 300),
        width=max(10*data.OTU.nunique(), 300),
        xrotation=90
    )

    kdims = ["OTU", "Group"]
    vdims = ["z_score"] + data.columns.drop(["z_score"]+kdims).tolist()

    hm = hv.HeatMap(
        data=data, kdims=kdims, vdims=vdims
    ).opts(**hm_opt).redim.values(Group=group_order, OTU=feature_order)

    hv.save(hm, f"clustermap_${otu_id}.html")
    """
}
