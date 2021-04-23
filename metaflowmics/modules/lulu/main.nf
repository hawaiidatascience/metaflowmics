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

    container "nakor/metaflowmics-script-env:0.0.1"
    // no conda repository for LULU

    input:
    tuple val(otu_id), path(matchlist), path(abundance), path(fasta)

    output:
    tuple val(otu_id), path("abundance_table-lulu-*.csv"), emit: csv
    tuple val(otu_id), path("*.fasta"), emit: fasta
    tuple val(otu_id), path("mapping_discarded*.txt"), emit: discarded, optional: true
    path "summary.csv", emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript

    library(lulu)
    library(data.table)
    library(seqinr)

    fast_table_load <- function(filename, drop=c(), row.label=1) {
        data <- as.data.frame(
          fread(filename, nThread=$task.cpus, drop=drop, header=T, blank.lines.skip=T)
        )
        if (row.label > 0) {
            rownames(data) <- data[[row.label]]
            data[, row.label] <- NULL
        }
        return(data)
    }

    if("${abundance.getExtension()}"=='shared') {
        otutab <- t(fast_table_load("$abundance", drop=c('label', 'numOtus'), row.label='Group'))
    } else {
        otutab <- fast_table_load("$abundance", row.label=1)        
    }

    matchList <- read.table("$matchlist", header=FALSE, col.names=c("OTU1", "OTU2", "pctIdentity"),
                            as.is=TRUE, check.names=F, stringsAsFactors=FALSE)
    fasta <- read.fasta("$fasta", seqtype="DNA", forceDNAtolower=F)

    if (dim(matchList)[1] > 0) {
        res <- lulu(as.data.frame(otutab), matchList, 
                    minimum_ratio_type=$params.lulu_min_ratio_type,
                    minimum_ratio=$params.lulu_min_ratio,
                    minimum_match=$params.lulu_min_match,
                    minimum_relative_cooccurence=$params.lulu_min_rel_cooccurence)
        fasta <- fasta[res\$curated_otus]

        otutab <- res\$curated_table
        write.table(res\$otu_map[res\$discarded_otus,], "mapping_discarded-${otu_id}.txt", quote=F)
    }

    write.csv(otutab, "abundance_table-lulu-${otu_id}.csv", quote=F)
    write.fasta(fasta, names(fasta), "sequences-${otu_id}.fasta")

    # Write counts
    summary <- cbind(
        rep('LULU-${otu_id}', ncol(otutab)),
        colnames(otutab),
        colSums(otutab),
        colSums(otutab > 0)
    )
    write.table(summary, "summary.csv", quote=F, sep=',', row.names=F, col.names=F)    

    writeLines(paste0(packageVersion('lulu')), "${software}.version.txt")
    """
}
