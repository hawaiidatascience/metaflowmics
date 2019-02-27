#!/usr/bin/env Rscript

library(ggplot2)
library(lulu)
library(dada2)
library(seqinr)
library(stringr)

extractChimeras <- function(derepFile,noChimeraFile,pairId)
{
    # extract non chimeric sequences from derep file
    derep <- readRDS(derepFile)
    seq.clean <- names(read.fasta(noChimeraFile,
                                  seqtype="DNA"))
    # get indices from no chimeric sequences
    indices <- as.numeric(sapply(seq.clean,
                                 function(x) str_extract(x,"[0-9]+")))

    # subset the indices from derep-class object
    derep[["uniques"]] <- derep[["uniques"]][indices]
    derep[["quals"]] <- derep[["quals"]][indices,]
    derep[["map"]] <- derep[["map"]][which(derep[["map"]] %in% indices)]
    
    # Handle special case where the longest sequence is removed
    seq.lengths <- sapply(names(derep[['uniques']]),nchar)
    derep[["quals"]] <- derep[["quals"]][,1:max(seq.lengths)]

    saveRDS(derep,paste0(pairId,"_nochimera.RDS"))
}

learnErrorRates <- function(noChimeraDerepFile,pairId)
{
    noChimeraDerep <- readRDS(noChimeraDerepFile)
    errors <- learnErrors(noChimeraDerep,
                          multithread=TRUE,
                          randomize=TRUE,
                          nbases=1e8)

    fig <- plotErrors(errors, nominalQ=TRUE)
    ggsave(paste0("R1_",pairId,".err.png"), plot=fig)

    saveRDS(errors, paste0(pairId,"_errors.RDS"))
}

dadaDenoise <- function(errorFile,derepFile,pairId)
{
    errors <- readRDS(errorFile)
    derep <- readRDS(derepFile)
    derep.denoised <- dada(derep, err=errors, multithread=TRUE)
    ## Extract abundance values from derep.denoised to name the fastas. Sample name also used
    seqIds <- paste0(pairId,";size=",as.numeric(getUniques(derep.denoised)))
    ## Write to fasta
    uniquesToFasta(derep.denoised,paste0(pairId,".dada.fasta"),seqIds)
    ## Save RDS object
    saveRDS(derep.denoised,paste0(pairId,".dada.RDS"))
}
    
esvTable <- function()
{
    ## Retrieve all the dada2 denoised files for merging
    denoisedFiles <- list.files(path=".", pattern="*.dada.RDS")
    pairIds <- as.character(sapply(denoisedFiles,
                                   function(x) unlist(strsplit(basename(x),"_R\\d{1}"))[1]))
    denoised <- lapply(denoisedFiles, function (x) readRDS(x))
    names(denoised) <- pairIds

    saveRDS(denoised,"denoised.RDS")
    
    seqtab <- makeSequenceTable(denoised)
    esv_ids <- paste0("esv_",c(1:dim(seqtab)[2]))
    uniquesToFasta(seqtab,"otus100_seq.fasta",ids=esv_ids)
    colnames(seqtab) <- esv_ids

    write.csv(t(seqtab),"otus100_table.csv",quote=F)

    ## Make summary file
    derepFiles <- lapply(list.files(path=".", pattern="*.derep.RDS"), function (x) readRDS(x))    
    chimeraFiles <- lapply(list.files(path=".", pattern="*.nochimera.RDS"), function (x) readRDS(x))

    summary.derep <- sapply(derepFiles, function(x) paste0(sum(x$uniques)," (",length(x$uniques)," uniques)"))        
    summary.chimera <- sapply(chimeraFiles, function(x) paste0(sum(x$uniques)," (",length(x$uniques)," uniques)"))    
    summary.denoised <- sapply(denoised, function(x) paste0(sum(x$denoised)," (",length(x$denoised)," uniques)"))

    summary <- data.frame("Sample" = pairIds,
                          "5-Dereplication" = summary.derep,
                          "6-ChimeraRemoval" = summary.chimera,
                          "7-Denoising" = summary.denoised,
                          check.names=F)
        
    write.table(summary, "count_summary.tsv", row.names=F, sep="\t")                    

}

luluCurate <- function(abundanceFile,matchListFile,threshold)
{
    # Load data
    otutab <- read.csv(abundanceFile,
                       header=TRUE,
                       as.is=TRUE,
		       check.names=F,
                       row.names=1)
    
    matchList <- read.table(matchListFile,
                            header=FALSE,
                            as.is=TRUE,
			    check.names=F,
                            col.names=c("OTU1","OTU2","pctIdentity"),
                            stringsAsFactors=FALSE)

    if (dim(matchList)[1] > 0) {
        ## Run Lulu
        curated <- lulu(otutab, matchList)

        ## Save result
        write.csv(curated$curated_table,
                  paste0("curated_table_",threshold,".csv"),
                  quote=F)
    
        write.table(curated$curated_otus,
                    paste0("curated_ids_",threshold,".csv"),
                    quote=F,
                    row.names=F,
                    col.names=F)
    } else {
        file.copy(abundanceFile, paste0("curated_table_",threshold,".csv"))
        write.table(rownames(otutab),
                    paste0("curated_ids_",threshold,".csv"),
                    quote=F,
                    row.names=F)
    }

}
