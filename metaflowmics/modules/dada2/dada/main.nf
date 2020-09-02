// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DADA2_DADA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "nakor/dada2:1.16"
    // conda (params.conda ? "bioconda::bioconductor-dada2=1.16 r-ggplot2" : null) // not working

    input:
    tuple val(meta), path(reads), path(errors)
    val options

    output:
    tuple val(meta), path("*-denoised.RDS"), emit: denoised
    tuple val(meta), path("*-derep.RDS"), emit: derep, optional: true
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(ggplot2)
    
    if (tools::file_ext("$reads") == 'RDS') {
        input <- readRDS("$reads")
    } else {
        input <- derepFastq("$reads")
        saveRDS(input, "${meta.id}-derep.RDS")
    }
    
    errors <- readRDS("$errors")
    denoised <- dada(input, err=errors, multithread=TRUE)
    saveRDS(denoised, "${meta.id}-denoised.RDS")

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")
    """
}
