// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process DADA2_ASSIGN_TAXONOMY {
    tag "${otu_id}"
    label "process_high"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h399db7b_1"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.18 conda-forge::r-ggplot2" : null)

    input:
    tuple val(otu_id), path(fasta)
    path db

    output:
    tuple val(otu_id), path("taxonomy*.csv"), emit: taxonomy
    tuple val(otu_id), path("confidence*.csv"), emit: confidence
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
    write.table(taxa, "taxonomy_${otu_id}.csv", quote=F, row.names=F, sep=",")

    scores <- cbind(seq_ids, assignments\$boot)
    write.table(scores, "confidence_${otu_id}.csv", quote=F, row.names=F, sep=",")

    writeLines(paste0(packageVersion("dada2")), "${software}.version.txt")
    """
}
