
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process HOLOVIEWS_SCATTER {
    tag "$meta.id"
    label "process_low"
    label "plot"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        meta:meta) }

    container "nakor/metaflowmics-python:0.0.2"
    conda (params.enable_conda ? "conda-forge::datatable pandas>=1" : null)

    input:
    tuple val(meta), file(mg)

    output:
    tuple val(meta), file("scatter*.html")

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
        data_r = data[data[rank].notnull()]

        if rank == "Class" and "Proteobacteria" in set(data_r.Phylum):
            data_r = data_r[data_r.Phylum == "Proteobacteria"]
    
        names = sorted(data_r[rank].unique())

        scatters = [hv.Scatter(
            data_r[data_r[rank] == name], "abundance", y_cols, label=name
        ).opts(**scatter_opt) for name in names]

        scatters = hv.Layout(scatters).cols(5)
        hv.save(scatters, f"scatterplot-{rank}_${meta.id}.html")
    """
}
