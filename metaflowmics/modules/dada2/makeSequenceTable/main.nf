// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DADA2_MAKESEQUENCETABLE {
    tag ""
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:"") }

    container "nakor/dada2:1.16"
    // conda (params.conda ? "bioconda::bioconductor-dada2=1.16 r-ggplot2" : null) // not working

    input:
    path merged_rds
    val options

    output:
    tuple val(100), path("dada2_ESVs-100.fasta"), emit: fasta
    tuple val(100), path("dada2_ESVs-100.csv"), emit: abundance
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env Rscript
    library(dada2)

    merged <- readRDS("${merged_rds}")
    merged.accepted <- lapply(merged, function(df) df[df\$accept,])
    esvTable <- makeSequenceTable(merged.accepted)

    print("Removing empty samples")
    not_empty_samples <- rowSums(esvTable)>0
    esvTable <- esvTable[not_empty_samples, ]
    
    write.csv(esvTable,"dada2_ESVs-100.csv")

    esv_ids <- sprintf("esv_%s", c(1:dim(esvTable)[2]))
    uniquesToFasta(esvTable, "dada2_ESVs-100.fasta", ids=esv_ids)
    colnames(esvTable) <- esv_ids

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")
    """
}
