#!/usr/bin/env Rscript

library(lulu)
library(dada2)
library(seqinr)
library(ShortRead)
library(stringr)

filterReads <- function(pairId,fwd,rev=NULL,
                        minLen=c(30,30),
                        maxEE=c(Inf,Inf),
                        truncLen=c(0,0),
                        rm.phix=TRUE,
                        truncQ=c(2,2)
                        )
{
    fwd.out <- paste(paste0(pairId,"1"),"trimmed.fastq",
                     sep="_")
    rev.out <- paste(paste0(pairId,"2"),"trimmed.fastq",
                     sep="_")

    if ( !is.null(rev) )
    {
        filterAndTrim(fwd, fwd.out, rev=rev, filt.rev=rev.out, compress=FALSE, truncQ=truncQ, minLen=minLen, maxEE=maxEE)
        saveRDS(readFastq(fwd.out)@id, paste0( pairId, "1.ids") )
        saveRDS(readFastq(rev.out)@id, paste0( pairId, "2.ids") )
    } else {
        filterAndTrim(fwd, fwd.out, compress=FALSE, truncQ=truncQ[1], minLen=minLen[1], maxEE=maxEE[1])
        saveRDS(readFastq(fwd.out)@id, paste0( pairId, "1.ids") )        
    }
            
}

learnErrorRates <- function(fastq,pairId)
{
    errors <- learnErrors(fastq,
                          multithread=TRUE,
                          randomize=TRUE,
                          nbases=1e8)
    pdf(paste0(pairId,".err.pdf"))
    plotErrors(errors, nominalQ=TRUE)
    dev.off()
    saveRDS(errors, paste0(pairId,"_errors.RDS"))
}

dadaDenoise <- function(errorFile,derepFile,pairId)
{
    errors <- readRDS(errorFile)
    derep <- derepFastq(derepFile)
    derep.denoised <- dada(derep, err=errors, multithread=TRUE, OMEGA_A=1e-80)
    # # Extract abundance values from derep.denoised to name the fastas. Sample name also used
    # seqIds <- paste0(pairId,";size=",as.numeric(getUniques(derep.denoised)))
    # Save RDS object
    saveRDS(derep,paste0(pairId,".derep.RDS"))
    saveRDS(derep.denoised,paste0(pairId,".dada.RDS"))
}

readsFromDenoised <- function(dadaOut, derep, pairId, readIds)
{
    denoised <- readRDS(dadaOut)
    sequences <- names(denoised$denoised)
    copies <- as.numeric(denoised$denoised)
    
    derepMap <- as.numeric(readRDS(derep)$map)
    denoisedMap <- as.numeric(denoised$map)
    
    reads <- do.call(function(i) { sequences[denoisedMap[i]] },
                     list(derepMap))

    ids <- readRDS(readIds)
    names(reads) <- unlist(strsplit(toString(ids),', '))

    writeFasta(reads,paste0(pairId,'.denoised.fasta'),width=1000)
}

    
esvTable <- function(pattern="*.dada.RDS")
{
    # Retrieve all the dada2 denoised files for merging
    mergerFiles <- list.files(path=".", pattern=pattern)
    pairIds <- as.character(sapply(mergerFiles,
                                   function(x) unlist(strsplit(basename(x),".dada"))[1]))
    mergers <- lapply(mergerFiles, function (x) readRDS(x))
    names(mergers) <- pairIds
    seqtab <- makeSequenceTable(mergers)
    esv_ids <- paste0("esv_",c(1:dim(seqtab)[2]))
    uniquesToFasta(seqtab,"esv_seq.fasta",ids=esv_ids)
    colnames(seqtab) <- esv_ids

    write.csv(t(seqtab),"esv_table.csv",quote=F)
}


luluCurate <- function(abundanceFile,matchListFile,threshold)
{
    otutab <- read.table(abundanceFile,
                         header=TRUE,
                         as.is=TRUE,
                         check.names=T,
                         row.names=2
                         )[-c(1,2)]
    otutab <- as.data.frame(t(otutab))
    
    matchList <- read.table(matchListFile,
                            header=FALSE,
                            as.is=TRUE,
                            stringsAsFactors=FALSE
                            )
    
    curated <- lulu(otutab, matchList, minimum_ratio_type="min", minimum_ratio=1, minimum_match=90, minimum_relative_cooccurence=0.95)
    
    write.csv(curated$curated_table,
              paste0("curated_table_",threshold,".csv"),
              quote=F)
    
    write.table(curated$curated_otus,
                paste0("curated_ids_",threshold,".csv"),
                quote=F,
                row.names=F,
                col.names=F)

}
