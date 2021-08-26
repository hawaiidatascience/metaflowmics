// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)


process BUILD_ASV_TABLE {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        meta:meta) }

    container "nakor/metaflowmics-r:0.0.1"    
    conda (params.enable_conda ? "conda-forge::r-stringr bioconda::bioconductor-dada2=1.18 conda-forge::r-seqinr" : null)

    input:
    path(rds)

    output:
    path "ASVs_duplicates_to_cluster.fasta", optional: true, emit: fasta_dup
    path "ASVs-100.fasta", emit: fasta    
    path "ASVs-100.{count_table,tsv}", emit: count_table
    // path "*.version.txt", emit: version

    script:
    // def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    library(dada2)
    library(seqinr)

    # Collect denoised reads
    denoised <- list.files(path=".", pattern="*.denoised.RDS")

    sample_names <- unname(sapply(denoised, function(x) gsub(".denoised.RDS", "", x)))
    merged <- lapply(denoised, function (x) readRDS(x))
    names(merged) <- sample_names
    
    # Retrieve merged object
    asv_table <- makeSequenceTable(merged)
    asv_ids <- sprintf("asv_%s", 1:dim(asv_table)[2])

    # Write ASV sequences
    uniquesToFasta(asv_table, "ASVs-100.fasta", ids=asv_ids)

    # Format count table
    count_table <- cbind(
        asv_ids, 
        colSums(asv_table),
        t(asv_table[rowSums(asv_table)>0,])
    )
    colnames(count_table) <- c("Representative_Sequence", "total", sample_names)

    if ("${params.format.toLowerCase()}" == "mothur") {
        # Write abundances
        write.table(count_table, file="ASVs-100.count_table", quote=F, sep="\\t",
                    row.names=F, col.names=T)
    } else {
        # Write abundances
        write.table(count_table[, -c(2)],"ASVs-100.tsv", quote=F, row.names=F, sep="\\t")

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
    """
}
