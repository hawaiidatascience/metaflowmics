// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process PHYLOSEQ_UNIFRAC {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process))}

    container "nakor/metaflowmics-r:0.0.2"
    conda (params.enable_conda ? "bioconda::bioconductor-phyloseq::1.34.0 conda-forge::r-data.table" : null)

    input:
    tuple val(meta), path(shared), path(tree)

    output:
    tuple val(meta), path("unifrac*.csv"), emit: csv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(phyloseq)

    abund <- data.frame(
        fread("$shared", drop=c("label", "numOtus"), sep="\\t"),
        row.names=1, check.names=F
    )

    tree <- read_tree("$tree")

    ps <- phyloseq(
        otu_table(as.matrix(abund), taxa_are_rows=FALSE),
        phy_tree(tree)
    )

    dists <- UniFrac(ps, weighted=("$params.unifrac"=="weighted"), parallel=TRUE)

    write.csv(as.matrix(dists), "unifrac_${params.unifrac}_${meta.id}.csv")

    writeLines(paste0(packageVersion("phyloseq")), "${software}.version.txt")
    """
}
