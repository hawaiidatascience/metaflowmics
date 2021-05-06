// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process SUMMARIZE_TABLE {
    tag "$step"
    label "process_low"

    container "nakor/metaflowmics-r:0.0.1"
    conda (params.enable_conda ? "conda-forge::r-data.table" : null)

    input:
    tuple val(step), val(otu_id), path(table)

    output:
    file('summary.csv')

    script:
    drop = 'NULL'
    sep = 'auto'
    rownames = 1
    taxa_are_rows = 'F'
    
    if (table.getExtension() == 'count_table') {
        drop = "'total'"
        sep = '\\t'
        rownames = 1
        taxa_are_rows = 'T'
    }
    """
    #!/usr/bin/env Rscript

    library(data.table)

    table <- data.frame(
        fread("$table", drop=$drop, sep='$sep'),
        row.names=$rownames, check.names=F
    )

    if (!$taxa_are_rows) {
        table <- t(table)
    }
    
    summary <- cbind(
        rep('$step', ncol(table)),
        rep('$otu_id', ncol(table)),
        colnames(table),
        colSums(table),
        colSums(table > 0)
    )

    write.table(summary, 'summary.csv', quote=F, row.names=F, col.names=F, sep=',')    
    """
    
}

process READ_TRACKING {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        meta:meta) }

    container "nakor/metaflowmics-r:0.0.1"
    conda (params.enable_conda ? "conda-forge::r-dplyr conda-forge::tidyr" : null)

    input:
    path counts

    output:
    file('summary-per-sample-per-step*.csv')

    script:
    """
    #!/usr/bin/env Rscript    

    library(stringr)
    library(dplyr)
    library(tidyr)

    data <- read.csv('summary.csv', header=F)
    colnames(data) <- c('step', 'otu_id', 'sample', 'total', 'nuniq')

    # Order the step according to total count and uniques
    col_order <- data %>% replace_na(list(nuniq=Inf)) %>%
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

    for (id in thresholds) {
        summary_i <- before_clustering %>%
          bind_rows(after_clustering %>% filter(otu_id==id)) %>%
          select(-otu_id) %>%
          pivot_wider(names_from=step, values_from=label)
        write.table(summary_i, sprintf('summary-per-sample-per-step-%s.csv', id), quote=F, row.names=F, sep=',')
    }
    """
}

process GET_SUBSAMPLING_THRESHOLD {
    label "process_low"

    container "nakor/metaflowmics-r:0.0.1"
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
    threshold <- floor(quantile(sample_sizes, ${params.subsampling_quantile}, names=F))

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
    container "nakor/metaflowmics-r:0.0.1"
    conda (params.enable_conda ? "conda-forge::r-data.table" : null)

    input:
    tuple val(otu_id), path(abundance), path(taxonomy)

    output:
    tuple val(otu_id), path("*.shared"), emit: shared
    tuple val(otu_id), path("*.taxonomy"), emit: taxonomy

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    abund <- data.frame(fread("$abundance"), row.names=1, check.names=F)

    if ($taxa_are_rows) {
        abund <- t(abund)
    }
    
    shared <- cbind(
        rep(1-$otu_id/100, nrow(abund)),
        rownames(abund),
        rep(ncol(abund), nrow(abund)),
        abund
    )
    colnames(shared) <- c('label', 'Group', 'numOtus', colnames(abund))
    write.table(shared, "${abundance.getName()}.shared", quote=F, sep='\\t', row.names=F)

    tax <- data.frame(fread("$taxonomy"), row.names=1, check.names=F)
    rank_names <- colnames(tax)

    tax <- cbind(
        rownames(tax), 
        rowSums(abund), 
        apply(tax, 1, function(x) paste(x, collapse=';'))
    )
    colnames(tax) <- c('OTU', 'Size', 'Taxonomy')

    write.table(tax, "${taxonomy.getName()}.taxonomy", quote=F, sep='\\t', row.names=F)
    """
}
