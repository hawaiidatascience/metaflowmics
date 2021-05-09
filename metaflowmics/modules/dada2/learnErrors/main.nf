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
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(stringr)
    library(ggplot2)

    inputs <- c("${reads.join('", "')}")
    
    for (input in inputs) {
        orient <- ifelse(str_detect(input, '_R2.'), '_R2',
                    ifelse(str_detect(input, '_R1.'), '_R1', ''))
        if (tools::file_ext(input) == 'RDS') {
            input <- readRDS(input)
        }

        tryCatch(
            expr={ errors <- learnErrors(input, multithread=TRUE, randomize=TRUE, nbases=1e7)
                   saveRDS(errors, sprintf("${meta.id}-errors%s.RDS", orient)) 
                   fig <- plotErrors(errors, nominalQ=TRUE)
                   ggsave(sprintf("${meta.id}-errorProfile%s.png", orient), plot=fig, type="cairo-png") },
            error=function(e) { print(e) }
        )
    }

    outputs <- list.files(pattern="-errors.*\\\\.RDS\$")
    if (length(outputs) < length(inputs)) {
        # Something went wrong (usually happens when there are too few reads
        # We discard the sample
        for (f in outputs) {
            file.remove(f)
        }
    }

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")    
    """
}
