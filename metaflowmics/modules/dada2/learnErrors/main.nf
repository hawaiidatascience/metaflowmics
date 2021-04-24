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
    tuple val(meta), path("*.RDS"), emit: rds
    path "*.png", emit: png
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    // def orient = (reads.getBaseName() =~/_(R[12])[\._]/)[0][1]
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(ggplot2)

    input <- "$reads"
    if (tools::file_ext("$reads") == 'RDS') {
        input <- readRDS("$reads")
    }

    errors <- learnErrors(input, multithread=TRUE, randomize=TRUE, nbases=1e8)
    saveRDS(errors, "${meta.id}-errors.RDS")

    fig <- plotErrors(errors, nominalQ=TRUE)
    ggsave("${meta.id}-errorProfile.png", plot=fig, type="cairo-png")

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")    
    """
}
