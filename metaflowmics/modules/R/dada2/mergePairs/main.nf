// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process DADA2_MERGEPAIRS {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "nakor/metaflowmics-r:0.0.2"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22 conda-forge::r-ggplot2 conda-forge::r-stringr conda-forge::r-seqinr conda-forge::r-dplyr conda-forge::r-tidyr" : null)

    input:
    path(derep)
    path(denoised)

    output:
    path "*.RDS", emit: rds
    tuple val(meta), path("ASVs.100.{count_table,tsv}"), emit: count_table
    tuple val(meta), path("ASVs.100.fasta"), emit: fasta
    tuple val(meta), path("ASVs_duplicates_to_cluster.fasta"), optional: true, emit: fasta_dup
    path "*_summary.tsv", emit: merge_summary, optional: true
    path "*.version.txt", emit: version

    script:
    meta = [id: "100"]
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript

    library(dada2)
    library(stringr)
    library(reshape2)
    library(dplyr)
    library(tidyr)
    library(seqinr)

    paired <- length(list.files(path=".", pattern="_R2.denoised.RDS")) > 1

    load_rds_dada2 <- function(type, orient="") {
        ext <- sprintf("%s.%s.RDS", orient, type)
        files <- list.files(path=".", pattern=ext)
        obj <- lapply(files, readRDS)
        names(obj) <- gsub(ext, "", files)
        return(obj)
    }

    ## Merge if paired_end (or merging was not done previously)
    if (paired) {
        merged <- mergePairs(
            derepF=load_rds_dada2("derep", "_R1"),
            dadaF=load_rds_dada2("denoised", "_R1"),
            derepR=load_rds_dada2("derep", "_R2"),
            dadaR=load_rds_dada2("denoised", "_R2"),
            minOverlap=${params.min_overlap},
            maxMismatch=${params.max_mismatch},
            returnReject=TRUE
        )
        saveRDS(merged, "merged.RDS")

        ## Mismatch summary
        mismatches <- lapply(merged, function(df) table(df\$nmismatch))
        summary <- reshape2::melt(mismatches, varnames="nmismatch") %>% 
          pivot_wider(values_from=value, names_from=L1) %>% 
          arrange(nmismatch)
        write.table(summary, "mismatch_summary.tsv", sep="\\t", quote=F, row.names=F)

        ## Removes the reads not matching the merging criteria
        merged <- lapply(merged, function(df) df[df\$accept,])
    } else { ## not paired
        merged <- load_rds_dada2("denoised", "")
        saveRDS(merged, "merged.RDS")
    }

    ## Make the ASV table
    asv_table <- makeSequenceTable(merged)
    asv_table <- asv_table[rowSums(asv_table)>0, ]
    sample_names <- rownames(asv_table)

    ## Write ASV sequences
    asv_ids <- sprintf("ASV_%s", c(1:dim(asv_table)[2]))
    uniquesToFasta(asv_table, "ASVs.100.fasta", ids=asv_ids)

    ## Compute count table
    count_table <- t(asv_table)
    rownames(count_table) <- asv_ids

    count_table <- cbind(asv_ids, rowSums(count_table), count_table)
    colnames(count_table) <- c("Representative_Sequence", "total", sample_names)

    if ("${params.format.toLowerCase()}" == "mothur") {
        write.table(count_table, file="ASVs.100.count_table", row.names=F, col.names=T, quote=F, sep="\\t")
    } else {
        write.table(count_table[, -c(2)], "ASVs.100.tsv", quote=F, row.names=F, sep="\\t")

        # Write duplicated fasta sequences with header formatted for VSEARCH
        list.fasta <- list()
        i = 1
        for(seq in colnames(asv_table)) {
            for(sample in rownames(asv_table)) {
                abd = asv_table[sample, seq]
                if(abd > 0) {
                    seq_id = sprintf("ASV_%s;sample=%s;size=%s", i, sample, abd)
                    list.fasta[seq_id] = seq
                }
            }
            i <- i+1
        }
        write.fasta(list.fasta, names=names(list.fasta), 
                    file.out='ASVs_duplicates_to_cluster.fasta')        
    }
    writeLines(paste0(packageVersion("dada2")), "${software}.version.txt")
    """
}
