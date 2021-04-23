// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process DADA2_DADA {
    tag "$meta.id"
    label 'process_low'
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

    for (read in c("${reads.join('"')}")) {
        orient <- ifelse(str_detect(read, '_R2.fastq'), '_R2', 
                    ifelse(str_detect(read, '_R1.fastq'), '_R1', ''))
        if (tools::file_ext(read) == 'RDS') {
            input <- readRDS(read)
        } else {
            input <- derepFastq(read)
            saveRDS(input, sprintf("${meta.id}-derep%s.RDS", orient))
        }
        errors <- readRDS("$errors")
        denoised <- dada(input, err=errors, multithread=TRUE)
        saveRDS(denoised, sprintf("${meta.id}-denoised%s.RDS", orient))
    }

    # Write counts
    counts <- getUniques(denoised)
    data <- sprintf("Denoising,${meta.id},%s,%s",sum(counts),sum(counts>0))
    write(data, "summary.csv")

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")
    """
}
