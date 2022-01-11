// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process LULU {
    tag "$meta.id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        pattern: "*.{csv,shared,fasta}",
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta) }

    container "nakor/metaflowmics-r:0.0.2"
    // no conda repository for LULU

    input:
    tuple val(meta), path(matchlist), path(abundance), path(fasta)

    output:
    tuple val(meta), path("abundance_table-lulu-*.{csv,shared}"), emit: abundance
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("mapping_discarded*.txt"), emit: discarded, optional: true
    path "*.summary", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ext = abundance.getExtension()
    """
    #!/usr/bin/env Rscript

    library(lulu)
    library(data.table)
    library(seqinr)

    fast_table_load <- function(filename, drop=c(), row.names=1) {
        data <- data.frame(
            fread(filename, nThread=$task.cpus, drop=drop, header=T),
            row.names=row.names, check.names=F
        )
        return(data)
    }

    if("$ext"=='shared') {
        otutab <- t(fast_table_load("$abundance", drop=c(1, 3), row.names=1))
    } else {
        otutab <- fast_table_load("$abundance", row.names=1)        
    }

    matchList <- read.table("$matchlist", header=FALSE, col.names=c("OTU1", "OTU2", "pctIdentity"),
                            as.is=TRUE, check.names=F, stringsAsFactors=FALSE)
    fasta <- read.fasta("$fasta", seqtype="DNA", forceDNAtolower=F)
    names(fasta) <- gsub('[\\t;].*', '', names(fasta))

    res <- lulu(as.data.frame(otutab), matchList, 
                minimum_ratio_type="$params.lulu_min_ratio_type",
                minimum_ratio=$params.lulu_min_ratio,
                minimum_match=$params.lulu_min_match,
                minimum_relative_cooccurence=$params.lulu_min_rel_cooccurence)
    otutab <- res\$curated_table
    fasta <- fasta[res\$curated_otus]

    if ("$ext"=='shared') {
        shared <- cbind(
            rep(1-${meta.otu_id}/100, ncol(otutab)),
            colnames(otutab),
            rep(nrow(otutab), ncol(otutab)),
            t(otutab)
        )
        colnames(shared) <- c('label', 'Group', 'numOtus', rownames(otutab))
        write.table(shared, "abundance_table-lulu-${meta.id}.shared", quote=F, sep='\\t', row.names=F)
    } else {
        abund <- cbind(
            rownames(otutab),
            otutab
        )
        colnames(abund) <- c("OTU", colnames(otutab))
        write.table(abund, "abundance_table-lulu-${meta.id}.csv", quote=F, sep=',', row.names=F)
    }

    write.table(res\$otu_map[res\$discarded_otus,], "mapping_discarded-${meta.id}.txt", quote=F)
    write.fasta(fasta, names(fasta), "sequences-${meta.id}.fasta")

    # Write counts
    summary <- cbind(
        rep("LULU", ncol(otutab)),
        rep("$meta.id", ncol(otutab)),
        colnames(otutab),
        colSums(otutab),
        colSums(otutab > 0)
    )
    write.table(summary, "lulu_${meta.id}.summary", quote=F, sep=',', row.names=F, col.names=F)    

    writeLines(paste0(packageVersion('lulu')), "${software}.version.txt")
    """
}
