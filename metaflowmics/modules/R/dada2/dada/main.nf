// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process DADA2_DADA {
    tag "$meta.id"
    label "${params.pool == 'F' ? 'process_medium' : 'process_high'}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22 conda-forge::r-stringr" : null)

    input:
    tuple val(meta), path(reads), path(errors)

    output:
    tuple val(meta), path("*.denoised.RDS"), emit: denoised
    tuple val(meta), path("*.derep.RDS"), emit: derep, optional: true
    path "summary.csv", optional: true, emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def suffix = meta.paired_end ? "_${meta.orient}" : ""
	def pool = params.pool == "pseudo" ? "'pseudo'" : params.pool
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(stringr)

    derep_files <- sort(list.files(".", pattern=".[^(errors)].RDS"))
    sample_names <- gsub("${suffix}\\\\.[a-z]+.RDS", "", derep_files)
    derep <- lapply(derep_files, readRDS)
    names(derep) <- sample_names

    err <- readRDS("$errors")

    denoised <- dada(derep, err=err, multithread=TRUE, pool=$pool)

    if (length(sample_names) == 1) {
        denoised <- list("$meta.id"=denoised)
    }

    sapply(names(denoised), function(x) saveRDS(denoised[[x]], sprintf("%s${suffix}.denoised.RDS", x)))

    # Write counts
    if ("$meta.orient" == "R1") {
        counts <- lapply(denoised, getUniques)
        abund <- sapply(counts, sum)
        richness <- sapply(counts, function(x) sum(x>0))
        data <- sprintf("denoising,,%s,%s,%s",sample_names,abund,richness)
        write(data, "summary.csv")
    }

    writeLines(paste0(packageVersion("dada2")), "${software}.version.txt")
    """
}
