// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DADA2_MERGEPAIRS {
    tag ""
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:"") }

    container "nakor/dada2:1.16"
    // conda (params.conda ? "bioconda::bioconductor-dada2=1.16 r-ggplot2" : null) // not working

    input:
    path(derep)
    path(denoised)
    val options

    output:
    path "reads_merged-dada2.RDS", emit: rds
    path "*.tsv", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env Rscript

    library(dada2)

    derepF <- list.files(path=".", pattern="*_R1-derep.RDS")
    derepR <- list.files(path=".", pattern="*_R2-derep.RDS")
    denoisedF <- list.files(path=".", pattern="*_R1-denoised.RDS")
    denoisedR <- list.files(path=".", pattern="*_R2-denoised.RDS")    
    
    inputs <- list(dadaF=lapply(denoisedF, readRDS),
                   derepF=lapply(derepF, readRDS),
                   dadaR=lapply(denoisedR, readRDS),
                   derepR=lapply(derepR, readRDS))

    params <- list(minOverlap=${params.min_overlap}, 
                   maxMismatch=${params.max_mismatch},
                   returnReject=TRUE)

    merged <- do.call(mergePairs,append(inputs, params))


    write_merge_summary <- function(merged, key, samples) {
        ## Get a summary of matches
        data <- lapply(merged, function(df) table(df[[key]]))
        indices <- sort(unique(unlist(lapply(data, function(x) as.numeric(names(x))))))

        summary <- as.data.frame(matrix(0, nrow=length(data), ncol=1+length(indices)))
        summary[, 1] <- samples
        colnames(summary) <- c('sample', indices)

        for (i in 1:nrow(summary)) {
            summary[i, names(data[[i]])] <- as.numeric(data[[i]])
        }
        write.table(summary, sprintf("%s_summary.tsv", key), sep='\\t', quote=FALSE, row.names=FALSE)
    }

    ## Get a summary of matches
    samples <- unname(sapply(denoisedF, function(x) gsub("_R1-denoised.RDS", "", x)))
    write_merge_summary(merged, 'nmatch', samples)
    write_merge_summary(merged, 'nmismatch', samples)

    ## Save data
    saveRDS(merged, "reads_merged-dada2.RDS")

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")    
    """
}
