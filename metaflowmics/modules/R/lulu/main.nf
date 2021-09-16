// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process LULU {
    tag "$otu_id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,
                                        publish_dir:getSoftwareName(task.process),
                                        meta:meta) }

    container "nakor/metaflowmics-r:0.0.1"
    // no conda repository for LULU

    input:
    tuple val(otu_id), path(matchlist), path(abundance), path(fasta)

    output:
    tuple val(otu_id), path("abundance_table-lulu-*.{csv,shared}"), emit: abundance
    tuple val(otu_id), path("*.fasta"), emit: fasta
    tuple val(otu_id), path("mapping_discarded*.txt"), emit: discarded, optional: true
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
            rep(1-${otu_id}/100, ncol(otutab)),
            colnames(otutab),
            rep(nrow(otutab), ncol(otutab)),
            t(otutab)
        )
        colnames(shared) <- c('label', 'Group', 'numOtus', rownames(otutab))
        write.table(shared, "abundance_table-lulu-${otu_id}.shared", quote=F, sep='\\t', row.names=F)
    } else {
        abund <- cbind(
            rownames(otutab),
            otutab
        )
        colnames(abund) <- c("OTU", colnames(otutab))
        write.table(abund, "abundance_table-lulu-${otu_id}.csv", quote=F, sep=',', row.names=F)
    }

    write.table(res\$otu_map[res\$discarded_otus,], "mapping_discarded-${otu_id}.txt", quote=F)
    write.fasta(fasta, names(fasta), "sequences-${otu_id}.fasta")

    # Write counts
    summary <- cbind(
        rep("LULU", ncol(otutab)),
        rep("$otu_id", ncol(otutab)),
        colnames(otutab),
        colSums(otutab),
        colSums(otutab > 0)
    )
    write.table(summary, "lulu_${otu_id}.summary", quote=F, sep=',', row.names=F, col.names=F)    

    writeLines(paste0(packageVersion('lulu')), "${software}.version.txt")
    """
}
