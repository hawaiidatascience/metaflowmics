// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process LULU {
    tag "$otu_id"
    label "process_high"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id: "$otu_id") }

    container "nakor/lulu"
    // no conda repository for LULU

    input:
    tuple val(otu_id), path(matchlist), path(abundance)
    val options

    output:
    tuple val(otu_id), path("OTU-ids-lulu-*.csv"), emit: ids
    tuple val(otu_id), path("abundance_table-lulu-*.csv"), emit: abundance
    tuple val(otu_id), path("mapping_discarded*.txt"), emit: discarded    
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    #!/usr/bin/env Rscript

    library(lulu)
    library(data.table)

    fast_table_load <- function(filename, cores=1, drop=c(), row.label=1) {
        data <- as.data.frame(data.table::fread(filename, nThread=cores, drop=drop, header=T, blank.lines.skip=T))
        if (row.label > 0) {
            rownames(data) <- data[[row.label]]
            data[, row.label] <- NULL
        }
        return(data)
    }

    if("${abundance.getExtension()}"=='shared'){
        otutab <- t(fast_table_load("$abundance", cores=$task.cpus, drop=c('label', 'numOtus'), row.label='Group'))
    } else {
        otutab <- fast_table_load("$abundance", cores=$task.cpus, row.label=1)        
    }

    matchList <- read.table("$matchlist", header=FALSE, col.names=c("OTU1", "OTU2", "pctIdentity"),
                            as.is=TRUE, check.names=F, stringsAsFactors=FALSE)

    if (dim(matchList)[1] > 0) {
        curated <- lulu(as.data.frame(otutab), matchList, 
                        minimum_ratio_type=$params.lulu_min_ratio_type,
                        minimum_ratio=$params.lulu_min_ratio,
                        minimum_match=$params.lulu_min_match,
                        minimum_relative_cooccurence=$params.lulu_min_rel_cooccurence)
    
    write.csv(curated\$curated_table, "abundance_table-lulu-${otu_id}.csv", quote=F)
    write.table(curated\$curated_otus, "OTU-ids-lulu-${otu_id}.csv", quote=F, row.names=F, col.names=F)
    write.table(curated\$otu_map[curated\$discarded_otus,], "mapping_discarded-${otu_id}.txt", quote=F)

    } else {
    write.csv(otutab, abundance_table_name)
    write.table(rownames(otutab), "OTU-ids-lulu-${otu_id}.csv", quote=F, row.names=F)

    empty.df <- read.table(text="", col.names=c('otu', 'parent_id'))
    write.table(empty.df, "mapping_discarded-${otu_id}.txt")
    }
 
    writeLines(paste0(packageVersion('lulu')), "${software}.version.txt")    
    """
}
