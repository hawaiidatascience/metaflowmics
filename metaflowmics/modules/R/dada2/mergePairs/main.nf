// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process DADA2_MERGEPAIRS {
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h399db7b_1"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.18 conda-forge::r-ggplot2" : null)

    input:
    path(derep)
    path(denoised)

    output:
    path "*.RDS", emit: rds
    path "ASVs-100.count_table", emit: count_table
    path "ASVs-100.fasta", emit: fasta
    path "*.tsv", emit: merge_summary, optional: true
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript

    library(dada2)

    derepF <- list.files(path=".", pattern="*-derep_R1.RDS")
    derepR <- list.files(path=".", pattern="*-derep_R2.RDS")
    denoisedF <- list.files(path=".", pattern="*-denoised_R1.RDS")
    denoisedR <- list.files(path=".", pattern="*-denoised_R2.RDS")

    sample_names <- unname(sapply(denoisedF, function(x) gsub("-denoised_R1.RDS", "", x)))

    merged <- mergePairs(
        derepF=lapply(derepF, readRDS),
        dadaF=lapply(denoisedF, readRDS),
        derepR=lapply(derepR, readRDS),
        dadaR=lapply(denoisedR, readRDS),
        minOverlap=${params.min_overlap},
        maxMismatch=${params.max_mismatch},
        returnReject=TRUE
    )
    names(merged) <- sample_names
    saveRDS(merged, "merged.RDS")

    # Merging summary
    write_merge_summary <- function(merged, key, samples) {
        ## Get a summary of matches
        data <- lapply(merged, function(df) table(df[[key]]))
        indices <- sort(unique(unlist(lapply(data, function(x) as.numeric(names(x))))))

        summary <- as.data.frame(matrix(0, nrow=length(data), ncol=1+length(indices)))
        summary[, 1] <- samples
        colnames(summary) <- c("sample", indices)

        for (i in 1:nrow(summary)) {
            summary[i, names(data[[i]])] <- as.numeric(data[[i]])
        }

        write.table(summary, sprintf("%s_summary.tsv", key), sep="\\t", quote=FALSE, row.names=FALSE)
    }

    write_merge_summary(merged, "nmatch", sample_names)
    write_merge_summary(merged, "nmismatch", sample_names)

    ## Removes the reads not matching the merging criteria
    merged <- lapply(merged, function(df) df[df\$accept,])

    ## Make the ASV table
    asv_table <- makeSequenceTable(merged)

    sample_names <- sample_names[rowSums(asv_table)>0]
    asv_table <- asv_table[sample_names, ]

    ## Write ASV sequences
    asv_ids <- sprintf("ASV_%s", c(1:dim(asv_table)[2]))
    uniquesToFasta(asv_table, "ASVs-100.fasta", ids=asv_ids)
    colnames(asv_table) <- asv_ids

    count_table <- cbind(asv_ids, colSums(asv_table), t(asv_table) )
    colnames(count_table) <- c("Representative_Sequence","total",sample_names)
    write.table(count_table, file="ASVs-100.count_table", row.names=F, col.names=T, quote=F, sep="\\t")
    writeLines(paste0(packageVersion("dada2")), "${software}.version.txt")
    """
}
