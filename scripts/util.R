#!/usr/bin/env Rscript

library(lulu)
library(dada2)
library(seqinr)
library(ShortRead)
library(stringr)

filterReads <- function(pairId,fwd,rev=NULL,
                        minLen=c(30,30),
                        maxEE=c(Inf,Inf),
                        truncLen=c(220,190),
                        rm.phix=TRUE,
                        truncQ=c(2,2)
                        )
{
    fwd.out <- paste(paste0(pairId,"_R1"),"trimmed.fastq",
                     sep="_")
    rev.out <- paste(paste0(pairId,"_R2"),"trimmed.fastq",
                     sep="_")

    if ( !is.null(rev) )
    {
        print(paste("Filtering paired-end reads with parameters",minLen,maxEE,truncLen,rm.phix,truncQ))
        filterAndTrim(fwd, fwd.out,
                      rev=rev, filt.rev=rev.out,
                      compress=FALSE,
                      truncLen=truncLen, truncQ=truncQ, minLen=minLen, maxEE=maxEE
                      )
        saveRDS(readFastq(fwd.out)@id, paste0( pairId, "_R1.ids") )
        saveRDS(readFastq(rev.out)@id, paste0( pairId, "_R2.ids") )

        pdf(paste0("qualityProfile_",pairId,".pdf"))
        plotQualityProfile(c(fwd,rev,fwd.out,rev.out))
        dev.off()
        
    } else {
        filterAndTrim(fwd, fwd.out,
                      compress=FALSE,
                      truncQ=truncQ[1], minLen=minLen[1], maxEE=maxEE[1], truncLen=truncLen
                      )
        saveRDS(readFastq(fwd.out)@id, paste0( pairId, "_R1.ids") )
        pdf(paste0("qualityProfile_",pairId,".pdf"))
        plotQualityProfile(c(fwd,fwd.out))
        dev.off()
    }
            
}

learnErrorRates <- function(fastq,pairId)
{
    errorsF <- learnErrors(fastq[1],
                           multithread=TRUE,
                           randomize=TRUE)
    saveRDS(errorsF, paste0(pairId,"_R1_errors.RDS"))
    
    if (length(fastq) > 1) {
        errorsR <- learnErrors(fastq[2],
                               multithread=TRUE,
                               randomize=TRUE)
        saveRDS(errorsR, paste0(pairId,"_R2_errors.RDS"))
    
        pdf(paste0("errorProfile_",pairId,".pdf"))
        plotErrors(c(errorsF,errorsR), nominalQ=TRUE)
        dev.off()
    } else {
        pdf(paste0("errorProfile_",pairId,".pdf"))
        plotErrors(errors, nominalQ=TRUE)
        dev.off()
    }
}

dadaDenoise <- function(errorFile,derepFile,pairId)
{
    errors <- readRDS(errorFile)
    derep <- derepFastq(derepFile)
    derep.denoised <- dada(derep, err=errors, multithread=TRUE)
    # # Extract abundance values from derep.denoised to name the fastas. Sample name also used
    # seqIds <- paste0(pairId,";size=",as.numeric(getUniques(derep.denoised)))
    # Save RDS object
    saveRDS(derep,paste0(pairId,".derep.RDS"))
    saveRDS(derep.denoised,paste0(pairId,".dada.RDS"))
}
    
esvTable <- function(minOverlap, maxMismatch, revRead)
{
    sample.names <- as.character(sapply( list.files(path=".",pattern="*_R1.derep.RDS"), 
                                         function(x) unlist(strsplit(x,"_R1",fixed=T))[1] )
    )

    derepF <- lapply(list.files(path=".",pattern="*_R1.derep.RDS"),readRDS)
    dadaF <- lapply(list.files(path=".",pattern="*_R1.dada.RDS"),readRDS)

    if (revRead == 1) {
        
        derepR <- lapply(list.files(path=".",pattern="*_R2.derep.RDS"),readRDS)
        dadaR <- lapply(list.files(path=".",pattern="*_R2.dada.RDS"),readRDS)

        merged <- mergePairs( dadaF, derepF,
                             dadaR, derepR,
                             minOverlap=minOverlap, maxMismatch=maxMismatch)
        saveRDS(merged,"dada_merged.RDS")

        esvTable <- makeSequenceTable(merged)
        esv.names <- paste0("sq",1:dim(esvTable)[2])

        uniquesToFasta(esvTable,"all.esv.fasta", ids=esv.names)

        esvTable <- cbind(esv.names, colSums(esvTable), t(esvTable) )

        colnames(esvTable) <- c("Representative_Sequence","total",sample.names)

        write.table(esvTable, file="all.esv.count_table", row.names=F, col.names=T, quote=F) 
                

    } else {
        print("Functionality not tested yet.")
        fasta_seq <- lapply( dadaF, function(x) names(x$denoised) )
        fasta_seq_uniq <- unique(cbind(unlist(fasta_seq)))
        esv.names <- paste0("sq",1:length(fasta_seq_uniq))

        write.fasta(fasta_seq_uniq,esv.names,"all.esv.fasta")
        
        esvTable <- as.data.frame(lapply(dadaF, function(x) x$denoised[fasta_seq_uniq]))

        esvTable[is.na(esvTable)] <- 0

        colnames(esvTable) <- sample.names

        esvTable <- cbind(esv.names, rowSums(esvTable), esvTable )
        colnames(esvTable) <- c("Representative_Sequence", "total", sample.names)
        
        write.table(esvTable, file="all.esv.count_table", row.names=F, col.names=T, quote=F)
    }
}

luluCurate <- function(abundanceFile,matchListFile,threshold)
{
    otutab <- read.table(abundanceFile,
                         header=TRUE,
                         as.is=TRUE,
                         check.names=F,
                         row.names=2
                         )[-c(1,2)]
    otutab <- as.data.frame(t(otutab))
    
    matchList <- read.table(matchListFile,
                            header=FALSE,
                            as.is=TRUE,
                            stringsAsFactors=FALSE
                            )
    
    curated <- lulu(otutab, matchList,
                    minimum_ratio_type="min",
                    minimum_ratio=1,
                    minimum_match=90,
                    minimum_relative_cooccurence=0.95)
    
    write.csv(curated$curated_table,
              paste0("curated_table_",threshold,".csv"),
              quote=F)
    
    write.table(curated$curated_otus,
                paste0("curated_ids_",threshold,".csv"),
                quote=F,
                row.names=F,
                col.names=F)

}
