// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process DADA2_FILTERANDTRIM {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h399db7b_1"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.18 conda-forge::r-ggplot2" : null)

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq, optional: true
    path("*.csv"), emit: summary
    path "*.png", emit: png, optional: true
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def rmphix = params.keep_phix ? "TRUE" : "FALSE"
    def n = meta.paired_end ? 2 : 1
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(ggplot2)

    filt_params <- data.frame(
        minLen=c($params.min_read_len),
        truncLen=c($params.trunc_len),
        truncQ=c($params.trunc_quality),
        maxEE=c($params.max_expected_error)
    )[1:$n,]

    if ("$meta.paired_end" == "true") {
        io <- list(fwd="${reads[0]}", filt="${meta.id}_R1.trimmed.fastq.gz",
                   rev="${reads[1]}", filt.rev="${meta.id}_R2.trimmed.fastq.gz")
    } else {
        io <- list(fwd="${reads[0]}", filt="${meta.id}.trimmed.fastq.gz")
    }

    params <- c(io, filt_params, list(rm.phix=${rmphix}))

    read_count <- do.call(filterAndTrim, params)[, "reads.out"]
    write(sprintf("qc,,$meta.id,%s,", read_count), "summary.csv")

    # Plot if we kept all reads
    if (file.exists(io[["filt"]])) {
        fig <- plotQualityProfile(io)
        ggsave("quality-profile_${meta.id}.png", plot=fig, type="cairo-png" )
    }

    writeLines(paste0(packageVersion("dada2")), "${software}.version.txt")
    """
}
