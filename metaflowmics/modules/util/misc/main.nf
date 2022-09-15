// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from "./functions"

options = initOptions(params.options)

process SUMMARIZE_TABLE {
    tag "$step"
    label "process_low"

    container "nakor/metaflowmics-r:0.0.2"
    conda (params.enable_conda ? "conda-forge::r-data.table" : null)

    input:
    tuple val(meta), val(step), path(table)

    output:
    tuple val(meta), file("summary.csv")

    script:
    drop = "NULL"
    sep = "auto"
    rownames = 1
    taxa_are_rows = "T"
	otu_id = meta.otu_id ?: ""
    
    if (table.getExtension() == "count_table") {
        drop = "'total'"
        sep = "\\t"
        rownames = 1
        taxa_are_rows = "T"
    }
    """
    #!/usr/bin/env Rscript

    library(data.table)

    table <- data.frame(
        fread("$table", drop=$drop, sep="$sep"),
        row.names=$rownames, check.names=F
    )

    if (!$taxa_are_rows) {
        table <- t(table)
    }
    
    summary <- cbind(
        rep("$step", ncol(table)),
        rep("$otu_id", ncol(table)),
        colnames(table),
        colSums(table),
        colSums(table > 0)
    )

    write.table(summary, "summary.csv", quote=F, row.names=F, col.names=F, sep=",")    
    """
    
}

process READ_TRACKING {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        meta:meta) }

    container "nakor/metaflowmics-r:0.0.2"
    conda (params.enable_conda ? "conda-forge::r-dplyr'=1.0.9' conda-forge::r-tidyr'=1.2' conda-forge::r-stringr" : null)

    input:
    path counts

    output:
    file("summary-per-sample-per-step*.csv")

    script:
	def outprefix = "summary-per-sample-per-step"
    """
    #!/usr/bin/env Rscript

    library(stringr)
    library(dplyr)
    library(tidyr)

    cols <- c("step", "otu_id", "sample", "total", "nuniq")
    data <- read.csv("$counts", header=F, col.names=cols)
    
    # Since mothur replaces "-" with "_", make sample names consistent
    data\$sample <- gsub("-", "_", data\$sample)

    # Order the step according to total count and uniques
    col_order <- data %>%
        mutate(nuniq=as.numeric(nuniq)) %>%
        replace_na(list(nuniq=Inf)) %>%
        group_by(step) %>% summarise(m1=sum(total), m2=sum(nuniq)) %>%
        arrange(desc(m1), desc(m2)) %>%
        pull(step) %>% as.character

    # Reshape the table into wide format
    summary <- data %>%
      mutate(
        step=factor(step, col_order),
        label=ifelse(is.na(nuniq), total, sprintf("%s (%s uniques)", total, nuniq))
      ) %>% 
      select(step, sample, label, otu_id) %>%
      arrange(step)

    # Split summary per clustering ids
    before_clustering <- summary %>% filter(is.na(otu_id))
    after_clustering <- summary %>% filter(!is.na(otu_id))
    thresholds <- after_clustering %>% pull(otu_id) %>% unique

    if (length(thresholds) > 0) {
        for (id in thresholds) {
            summary_i <- before_clustering %>%
              bind_rows(after_clustering %>% filter(otu_id==id)) %>%
              select(-otu_id) %>%
              pivot_wider(names_from=step, values_from=label, values_fill="0")
            write.table(summary_i, sprintf("${outprefix}.%s.csv", id), quote=F, row.names=F, sep=",")
        }
    } else { # reads were all filtered before clustering
        summary <- summary %>% 
            select(-otu_id) %>% 
            pivot_wider(names_from=step, values_from=label, values_fill="0")
        write.table(summary, "${outprefix}.csv", quote=F, row.names=F, sep=",")
    }
    """
}

process GET_SUBSAMPLING_THRESHOLD {
    label "process_low"
	tag "$params.subsampling_quantile"

    container "nakor/metaflowmics-r:0.0.2"
    conda (params.enable_conda ? "conda-forge::r-data.table" : null)

    input:
    path count

    output:
    stdout

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    table <- data.frame(fread("$count", drop=c(2)), row.names=1, check.names=F)
    sample_sizes <- colSums(table)
    threshold <- floor(quantile(sample_sizes, $params.subsampling_quantile, names=F))

    if (threshold < $params.min_subsampling) {
        threshold <- ifelse($params.min_subsampling < max(sample_sizes), 
                            $params.min_subsampling,
                            threshold)
    }

    cat(threshold)
    """
}

process CONVERT_TO_MOTHUR_FORMAT {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }	

    container "nakor/metaflowmics-r:0.0.2"
    conda (params.enable_conda ? "conda-forge::r-data.table" : null)

    input:
    tuple val(meta), path("abundance.tsv"), path("taxonomy.csv")

    output:
    tuple val(meta), path("*.shared"), emit: shared
    tuple val(meta), path("*.taxonomy"), emit: taxonomy

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    abund <- data.frame(fread("abundance.tsv"), row.names=1, check.names=F)

    if ($params.taxa_are_rows) {
        abund <- t(abund)
    }
    
    shared <- cbind(
        rep(1-$meta.otu_id/100, nrow(abund)),
        rownames(abund),
        rep(ncol(abund), nrow(abund)),
        abund
    )
    colnames(shared) <- c("label", "Group", "numOtus", colnames(abund))
    write.table(shared, "OTUs.${meta.id}.shared", quote=F, sep="\\t", row.names=F)

    tax <- data.frame(fread("taxonomy.csv"), row.names=1, check.names=F)
    rownames(tax) <- gsub(";.*", "", rownames(tax))
    rank_names <- colnames(tax)

    tax <- cbind(
        rownames(tax), 
        colSums(abund), 
        apply(tax, 1, function(x) paste(x, collapse=";"))
    )
    colnames(tax) <- c("OTU", "Size", "Taxonomy")
    tax[, "Taxonomy"] <- gsub("[a-z]__", "", tax[, "Taxonomy"])

    write.table(tax, "OTUs.${meta.id}.taxonomy", quote=F, sep="\\t", row.names=F)
    """
}
