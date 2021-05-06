
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process HOLOVIEWS_SCATTER {
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
    tuple val(otu_id), file("scatter*.html")

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env python

    import pandas as pd
    import holoviews as hv

    hv.extension("bokeh")

    scatter_opt = dict(
        size=8, alpha=0.8, 
        tools=["hover"], 
        axiswise=True
    )

    data = pd.read_csv("${mg}")
    numeric_data = data.groupby("OTU")["count"].agg(prevalence=lambda x: sum(x>0), abundance=sum)
    metadata = data.groupby("OTU").first().drop(columns=["Group", "count", "relabund"])
    data = pd.concat([numeric_data, metadata], axis=1)

    y_cols = ["prevalence"] + metadata.columns.tolist()

    for rank in ["Phylum", "Class"]:
        if rank == "Class":
            data = data[data.Phylum == "Proteobacteria"]
    
        names = sorted(data[rank].unique())

        scatters = [hv.Scatter(
            data[data[rank] == name], "abundance", y_cols, label=name
        ).opts(**scatter_opt) for name in names]

        scatters = hv.Layout(scatters).cols(5)
        hv.save(scatters, f"scatterplot-{rank}_${otu_id}.html")
    """
}
