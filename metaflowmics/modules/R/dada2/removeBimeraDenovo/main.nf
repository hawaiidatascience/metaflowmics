// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process DADA2_CHIMERA {
    tag "$meta.id"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta, publish_by_meta:["id"]) }

    container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h399db7b_1"
    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.18 conda-forge::r-ggplot2" : null)

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*.RDS"), emit: rds    
    path("*.csv"), emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript

    library(dada2)

    derep <- readRDS("$contigs")
    chimera <- isBimeraDenovo(derep, verbose=T)

    kept <- (1:length(chimera))[!chimera]

    # Remove uniques sequences identified as chimeric
    derep[["uniques"]] <- derep[["uniques"]][kept]

    # Find new maximum read length (in case we removed the longest)
    max_len <- max(sapply(names(derep[["uniques"]]), nchar))
    # Subset the qualities
    derep[["quals"]] <- derep[["quals"]][kept, 1:max_len]

    # Subset the read -> unique mapping
    derep[["map"]] <- derep[["map"]][which(derep[["map"]] %in% kept)]

    # Write counts
    counts <- getUniques(derep)
    data <- sprintf("chimera,,${meta.id},%s,%s",sum(counts),sum(counts>0))
    write(data, "summary.csv")

    saveRDS(derep, "${meta.id}.nochim.RDS")

    writeLines(paste0(packageVersion("dada2")), "${software}.version.txt")
    """
}
