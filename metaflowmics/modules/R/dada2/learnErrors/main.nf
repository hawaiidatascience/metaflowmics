// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process DADA2_LEARNERRORS {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:['id']) }

    container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h399db7b_1"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.18 conda-forge::r-ggplot2" : null)

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.RDS"), emit: rds, optional: true
    path "*.png", emit: png, optional: true
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def suffix = meta.paired_end ? "_${meta.orient}" : ""
    def outprefix = "${meta.id}${suffix}"
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(stringr)
    library(ggplot2)
    
    files <- sort(list.files(".", pattern="${suffix}.*.RDS"))
    sample_names <- gsub("${suffix}.*.RDS", "", files)
    reads <- lapply(files, readRDS)
    names(reads) <- sample_names

    tryCatch(
        expr={ errors <- learnErrors(reads, multithread=TRUE, randomize=TRUE, nbases=1e7)
               saveRDS(errors, "${outprefix}.errors.RDS") 
               fig <- plotErrors(errors, nominalQ=TRUE)
               ggsave("profile_${outprefix}.png", plot=fig, type="cairo-png") },
        error=function(e) { print(e) }            
    )

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")    
    """
}
