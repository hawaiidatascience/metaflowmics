// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process DADA2_ASSIGN_TAXONOMY {
    tag "${meta.id}"
    label "process_high"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22 conda-forge::r-ggplot2" : null)

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("taxonomy*.csv"), emit: taxonomy
    tuple val(meta), path("confidence*.csv"), emit: confidence
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    library(dada2)
    library(ShortRead)

    sr <- readFasta("$fasta")
    seq_ids <- as.character(id(sr))
    seq <- as.character(sread(sr))
    assignments <- assignTaxonomy(
        setNames(seq, seq_ids), "$db",
        minBoot=$params.tax_confidence,
        tryRC=TRUE,
        outputBootstraps=TRUE
    )
    taxa = cbind(seq_ids, assignments\$tax)
    write.table(taxa, "taxonomy_${meta.id}.csv", quote=F, row.names=F, sep=",")

    scores <- cbind(seq_ids, assignments\$boot)
    write.table(scores, "confidence_${meta.id}.csv", quote=F, row.names=F, sep=",")

    writeLines(paste0(packageVersion("dada2")), "${software}.version.txt")
    """
}
