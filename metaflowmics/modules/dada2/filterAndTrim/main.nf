// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DADA2_FILTERANDTRIM {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "nakor/dada2:1.16"
    // conda (params.conda ? "bioconda::bioconductor-dada2=1.16 r-ggplot2" : null) // not working

    input:
    tuple val(meta), path(reads)
    val options

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "*.png", emit: png
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)

    def rmphix = params.keepPhix ? "TRUE" : "FALSE"
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(ggplot2)

    params <- list(
        minLen=c($params.minLen), truncLen=c($params.truncLen), truncQ=c($params.truncQ), 
        maxEE=c($params.maxEE), rm.phix=${rmphix}, compress=TRUE
    )

    io <- list(fwd="${reads[0]}", filt="${meta.id}_R1-trimmed.fastq.gz")

    if ("${meta.paired}" == "true") {
        io[["rev"]] = "${reads[1]}"
        io[["filt.rev"]] = "${meta.id}_R2-trimmed.fastq.gz"
    }

    do.call(filterAndTrim, append(io, params))

    fig <- plotQualityProfile(io)
    ggsave("quality-profile_${meta.id}.png", plot=fig, type="cairo-png" )

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")    
    """
}
