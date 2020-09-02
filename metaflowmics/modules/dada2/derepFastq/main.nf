// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DADA2_DEREPFASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "nakor/dada2:1.16"
    // conda (params.conda ? "bioconda::bioconductor-dada2=1.16 r-ggplot2" : null) // not working

    input:
    tuple val(meta), path(reads)
    val options

    output:
    tuple val(meta), path("*.RDS"), emit: rds
    tuple val(meta), path("*.fasta"), emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env Rscript

    library(dada2)

    derep <- derepFastq("$reads")
    saveRDS(derep,"${meta.id}-derep.RDS")
	uniquesToFasta(derep,"${meta.id}-derep.fasta")

    writeLines(paste0(packageVersion('dada2')), "${software}.version.txt")    
    """
}
