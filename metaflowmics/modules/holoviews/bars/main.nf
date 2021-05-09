
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process HOLOVIEWS_BARS {
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
    tuple val(otu_id), file("barplot*.html")

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env python

    import pandas as pd
    import holoviews as hv

    hv.extension("bokeh")

    data = pd.read_csv("${mg}").drop(columns=["OTU"])

    n_samples = data.Group.nunique()
    
    legend_opts = dict(
        label_text_font_size="12px",
        glyph_height=12,
        spacing=0,
        label_height=12
     )

    plot_opts = dict(
        width=800, height=max(20*n_samples, 300),
        stacked=True, tools=["hover"], 
        invert_axes=True,
        legend_position="right",
        legend_offset=(30, 0),
        legend_opts=legend_opts
    )

    filler = f"Others(<{${params.min_abund*100}:g}%)"

    def combine(group, filler=filler):
        out = pd.Series(dict(
            (col, vals.sum()) if col in {"relabund", "count"} 
            else (col, vals.iloc[0])
            for (col, vals) in group.iteritems()
        ))

        if out["relabund"] < $params.min_abund:
            tax_cols = group.drop(columns=["Group", "relabund", "count"]).columns
            out[tax_cols] = filler

        return out

    for rank in ["Phylum", "Class"]:

        tax_cols = data.loc[:, "Kingdom":rank].columns.tolist()

        df = (data.loc[:, "Group":rank]
            .groupby(["Group", rank]).apply(combine)
            .reset_index(drop=True))

        if rank == "Class" and "Proteobacteria" in set(data.Phylum):
            df = df[df.Phylum == "Proteobacteria"]
            df.relabund = df.groupby("Group").relabund.transform(lambda x: x/sum(x))

        sorted_taxa = df.groupby(rank).relabund.sum().sort_values().index[::-1]
        sorted_groups = sorted(df.Group.unique())[::-1]

        if filler in sorted_taxa:
            sorted_taxa = sorted_taxa.drop(filler).append(pd.Index([filler]))

        bars = (
            hv.Bars(df.reset_index(), kdims=["Group", rank],
                    vdims=["relabund", "count"] + tax_cols[:-1])
            .redim.values(**{rank: sorted_taxa, "Group": sorted_groups})
            .opts(**plot_opts)
        )
        hv.save(bars, f"barplot-{rank}_${otu_id}.html")

    """
}
