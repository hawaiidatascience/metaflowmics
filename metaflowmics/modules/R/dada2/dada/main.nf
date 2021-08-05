// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process DADA2_DADA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:['id']) }

    container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h399db7b_1"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.18 conda-forge::r-stringr" : null)

    input:
    tuple val(meta), path(reads), path(errors)

    output:
    tuple val(meta), path("*-denoised*.RDS"), emit: denoised
    tuple val(meta), path("*-derep*.RDS"), emit: derep, optional: true
    path "summary.csv", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(stringr)

    read_files <- c("${reads.join('", "')}")
    error_files <- c("${errors.join('", "')}")

    for (i in 1:length(read_files)) {
        rd_file <- read_files[i]
        err_file <- error_files[i]

        orient <- ifelse(str_detect(err_file, '_R2.RDS'), '_R2', 
                    ifelse(str_detect(err_file, '_R1.RDS'), '_R1', ''))
        if (tools::file_ext(rd_file) == 'RDS') {
            input <- readRDS(rd_file)
        } else {
            input <- derepFastq(rd_file)
            saveRDS(input, sprintf("${meta.id}-derep%s.RDS", orient))
        }
        errors <- readRDS(err_file)
        denoised <- dada(input, err=errors, multithread=TRUE)
        saveRDS(denoised, sprintf("${meta.id}-denoised%s.RDS", orient))
    }

    # Write counts
    counts <- getUniques(denoised)
    data <- sprintf("denoising,,${meta.id},%s,%s",sum(counts),sum(counts>0))
    write(data, "summary.csv")

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")
    """
}
