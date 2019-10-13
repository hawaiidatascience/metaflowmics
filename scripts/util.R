#!/usr/bin/env Rscript

library(stringr)

library(ShortRead)
library(seqinr)

library(dada2)
library(lulu)

library(ggplot2)

filterReads <- function(pairId,fwd,rev=NULL,
                        minLen=c(30,30),
                        maxEE=c(Inf,Inf),
                        truncLen=c(220,190),
                        rm.phix=TRUE,
                        truncQ=c(2,2)
                        )
{
    fwd.out <- sprintf("%s_R1_trimmed.fastq.gz",pairId)
    rev.out <- sprintf("%s_R2_trimmed.fastq.gz",pairId)

    params <- list(truncQ=truncQ, truncLen=truncLen, minLen=minLen, maxEE=maxEE,rm.phix=rm.phix,compress=TRUE)
    
    if (is.null(rev)) {
        files.io <- list(fwd=fwd, filt=fwd.out)
    } else {
        files.io <- list(fwd=fwd, rev=rev, filt=fwd.out, filt.rev=rev.out)
    }

    print(append(files.io,params))
        
    ## Apply dada2's filterAndTrim
    do.call(filterAndTrim, append(files.io,params))

    if(!file.exists(files.io$filt)) {
        plotQualityProfile(files.io[c(1,2)])
        return()
    }
    ## Plot error profiles
    fig <- plotQualityProfile(files.io)
    ggsave( sprintf("qualityProfile_%s.png",pairId), plot=fig, type = "cairo-png" )
}

extractChimeras <- function(derepFile,noChimeraFile,pairId)
{
    ## Extract non chimeric sequences from derep file
    derep <- readRDS(derepFile)
    seq.clean <- names(read.fasta(noChimeraFile,
                                  seqtype="DNA"))
    ## Get indices from no chimeric sequences
    indices <- as.numeric(sapply(seq.clean,
                                 function(x) str_extract(x,"[0-9]+")))

    ## Subset the indices from derep-class object
    derep[["uniques"]] <- derep[["uniques"]][indices]
    derep[["quals"]] <- derep[["quals"]][indices,]
    derep[["map"]] <- derep[["map"]][which(derep[["map"]] %in% indices)]
    
    ## Handle special case where the longest sequence is removed
    seq.lengths <- sapply(names(derep[['uniques']]),nchar)
    derep[["quals"]] <- derep[["quals"]][,1:max(seq.lengths)]

    saveRDS(derep,paste0(pairId,"_R12_nochimera.RDS"))
}

learnErrorRates <- function(inputFiles,pairId,rds=FALSE)
{
    prefix.out = pairId
    if(length(inputFiles) > 1) {
        prefix.out <- sprintf("%s_%s",prefix.out,c('R1','R2'))
    }

    for( i in 1:length(inputFiles) ) {
        inputFile = inputFiles[i]
        if(rds) {
            inputFile <- readRDS(inputFile)
        }
        
        errors <- learnErrors(inputFile,multithread=TRUE,randomize=TRUE,nbases=1e8)
        saveRDS(errors, sprintf("%s_errors.RDS",prefix.out[i]))
        
        fig <- plotErrors(errors, nominalQ=TRUE)
        ggsave(sprintf("%s_errorProfile.png",prefix.out[i]), plot=fig, type="cairo-png")
    }
}

dadaDenoise <- function(errorFile,derepFile,pairId,derep.rds=FALSE,paired=TRUE)
{
    errors <- readRDS(errorFile)

    if(derep.rds) {
        derep <- readRDS(derepFile)
    } else {
        derep <- derepFastq(derepFile)
    }

    if(paired) {
        saveRDS(derep,sprintf("%s_derep.RDS",pairId))
    }
    
    denoised <- dada(derep, err=errors, multithread=TRUE)
    ## Save RDS object
    saveRDS(denoised,sprintf("%s_dada.RDS",pairId))
}

getCounts <- function(attribute) {
    getCounts_attr <- function(x) {
        #' Get the number of total and unique sequences
        #' on [attribute] attribute of object [x]
        total_count <- sum(x[[attribute]])
        uniq_count <- sum(x[[attribute]]>0)
        return(sprintf("%s (%s uniques)",total_count,uniq_count))
    }
    return(getCounts_attr)
}

mergeReads <- function(minOverlap, maxMismatch)
{
    derepF <- list.files(path=".",pattern="*_R1.derep.RDS")
    derepR <- list.files(path=".",pattern="*_R2.derep.RDS")
    denoisedF <- list.files(path=".",pattern="*_R1.dada.RDS")
    denoisedR <- list.files(path=".",pattern="*_R2.dada.RDS")    
    
    inputs <- list(dadaF=lapply(denoisedF, readRDS),
                   derepF=lapply(derepF, readRDS),
                   dadaR=lapply(denoisedR, readRDS),
                   derepR=lapply(derepR, readRDS))
    params <- list(minOverlap=minOverlap,maxMismatch=maxMismatch)

    merged <- do.call(mergePairs,append(inputs,params))

    saveRDS(merged,"reads_merged.RDS")

    return(merged)
}

makeSummary <- function(steps,patterns,attributes) {
    sample.names <- as.character(sapply(
        list.files(pattern=patterns[1]),
        function(x) unlist(strsplit(x,"_R1",fixed=T))[1]
    ))
    summary <- data.frame(matrix(ncol=1+length(steps), nrow=length(sample.names)))
    colnames(summary) <- c("Sample",steps)
    summary['Sample'] = sample.names
    
    for( i in 1:length(patterns) ) {
        files <- list.files(pattern=patterns[i])
        data <- lapply(files, readRDS)
        summary[steps[i]] = sapply(data, getCounts(attributes[i]))
    }

    return(summary)
}
    
esvTable <- function(minOverlap=20, maxMismatch=1, paired=TRUE)
{
    if(paired) {
        merged <- mergeReads(minOverlap, maxMismatch)
    } else {
        ## Retrieve all the dada2 denoised files for merging
        denoisedFiles <- list.files(pattern="*_dada.RDS")
        pairIds <- as.character(sapply(denoisedFiles,
                                       function(x) unlist(strsplit(basename(x),"_"))[1]))
        merged <- lapply(denoisedFiles, function (x) readRDS(x))
        names(merged) <- pairIds
    }

    esvTable <- makeSequenceTable(merged)
    write.csv(esvTable,"raw_sequence_table.csv")
    
    esv_ids <- sprintf("esv_%s",c(1:dim(esvTable)[2]))
    uniquesToFasta(esvTable,"all_esv.fasta",ids=esv_ids)
    colnames(esvTable) <- esv_ids

    ## Make summary file
    if( length(list.files(pattern="*_nochimera.RDS")) > 0 ) {
        ## ITS pipeline
        steps = c("4-Dereplication","5-ChimeraRemoval","6-Denoising")
        patterns = c("*_R1.*_derep.RDS","*_R1.*_nochimera.RDS","*_R1.*_dada.RDS")
        attributes = c("uniques","uniques","denoised")
        summary <- makeSummary(steps,patterns,attributes)

        write.csv(t(esvTable),"clustering_100.csv",quote=F)
    } else {
        ## 16S pipeline
        steps = c("3.1-Dereplication","3.2-Denoising")
        patterns = c("*_R1.*_derep.RDS","*_R1.*_dada.RDS")
        attributes = c("uniques","denoised")        
        summary <- makeSummary(steps,patterns,attributes)

        if(paired) {
            summary['4-ESV'] = sapply(merged, getCounts('abundance'))
        } else {
            summary['4-ESV'] = summary['3.2-Denoising']
        }
        
        ## write in .count_table format for Mothur
        esvTable <- cbind(esv_ids, colSums(esvTable), t(esvTable) )
        colnames(esvTable) <- c("Representative_Sequence","total",unlist(summary['Sample']))
        write.table(esvTable, file="all_esv.count_table", row.names=F, col.names=T, quote=F, sep="\t")
    }

    write.table(summary, "count_summary.tsv", row.names=F, sep="\t")                    
}

mergeFastaForVSEARCH <- function(filename)
{
    #' Retrieve the sequences from abundance table and create a fasta file
    #' The header of each sequence includes its sample and its abundance
    
    table <- read.csv(filename,row.names=1)
    list.fasta <- list()
    i = 1

    for(seq in colnames(table)) {
        in.sample = rownames(table)[which(table[,seq]>0, arr.ind=TRUE)]
        for(sample in in.sample) {
            seq.id = sprintf("esv_%s;sample=%s;size=%s",i,sample,table[sample,seq])
            list.fasta[seq.id] = seq
            i = i+1
        }
    }
    write.fasta(list.fasta, names=names(list.fasta), file.out='esv_merged_to_cluster.fasta')
}

luluCurate <- function(abundanceFile,matchListFile,threshold,min_ratio_type,min_ratio,min_match,min_rel_cooccurence)
{
    if(tools::file_ext(abundanceFile)=='shared'){
        otutab <- read.table(abundanceFile,
                             header=TRUE, row.names=2, colClasses=c("Group"="character"),
                             as.is=TRUE,check.names=F)[,-c(1,2)]
        otutab <- as.data.frame(t(otutab))
        abundance_table_name = sprintf("lulu_table_%s.csv",threshold)
    } else {
        otutab <- read.csv(abundanceFile,
                           header=TRUE, row.names=1,
                           as.is=TRUE, check.names=F)
        
        abundance_table_name = sprintf("abundance_table_%s.csv",threshold)
    }
        
    curated_ids_name = sprintf("lulu_ids_%s.csv",threshold)
    
    matchList <- read.table(matchListFile,
                            header=FALSE, col.names=c("OTU1","OTU2","pctIdentity"),
                            as.is=TRUE, check.names=F, stringsAsFactors=FALSE)

    if (dim(matchList)[1] > 0) {
        ## Run Lulu
        curated <- lulu(otutab, matchList,
                        minimum_ratio_type=min_ratio_type,
                        minimum_ratio=min_ratio,
                        minimum_match=min_match,
                        minimum_relative_cooccurence=min_rel_cooccurence)
    
        ## Save result
        write.csv(curated$curated_table,
                  abundance_table_name,
                  quote=F)
    
        write.table(curated$curated_otus,
                    curated_ids_name,
                    quote=F,
                    row.names=F,
                    col.names=F)
    } else {
        write.csv(otutab, abundance_table_name)
        write.table(rownames(otutab),
                    curated_ids_name,
                    quote=F,
                    row.names=F)
    }

}
