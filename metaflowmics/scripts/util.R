#!/usr/bin/env Rscript

library(stringr)
library(doParallel)

library(ShortRead)
library(seqinr)
library(dada2)
library(lulu)
library(ape)
library(phyloseq)

library(ggplot2)

fast_table_load <- function(filename, cores=1, drop=c(), row.label=1) {
    data <- as.data.frame(data.table::fread(filename, nThread=cores, drop=drop, header=T, blank.lines.skip=T))

    if (row.label > 0) {
        rownames(data) <- data[[row.label]]
        data[, row.label] <- NULL
    }
    return(data)
}

filter_reads <- function(pairId,fwd,rev=NULL,
                        minLen=30,
                        maxEE=c(Inf,Inf),
                        truncLen=c(220,190),
                        rm.phix=TRUE,
                        truncQ=2
                        )
{
    fwd.out <- sprintf("%s_R1_trimmed.fastq.gz",pairId)
    rev.out <- sprintf("%s_R2_trimmed.fastq.gz",pairId)

    params <- list(truncQ=truncQ, truncLen=truncLen, minLen=minLen, maxEE=maxEE,rm.phix=rm.phix,compress=TRUE)
    
    if (is.null(rev)) {
        files.io <- list(fwd=fwd, filt=fwd.out)

        for (name in c('minLen', 'maxEE', 'truncLen', 'truncQ')) {
            params[[name]] <- params[[name]][1]
        }
    } else {
        files.io <- list(fwd=fwd, rev=rev, filt=fwd.out, filt.rev=rev.out)
    }

    print(append(files.io,params))
        
    ## Apply dada2's filterAndTrim
    do.call(filterAndTrim, append(files.io,params))

    if (!file.exists(files.io$filt)) {
        if (is.null(rev)) {
            plotQualityProfile(files.io[[1]])
        } else {
            plotQualityProfile(files.io[c(1,2)])
        }
        return()
    }
    ## Plot error profiles
    fig <- plotQualityProfile(files.io)
    ggsave( sprintf("qualityProfile_%s.png",pairId), plot=fig, type = "cairo-png" )
}

extract_chimeras <- function(derepFile,noChimeraFile,pairId)
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

learn_error_rates <- function(inputFiles,pairId,rds=FALSE)
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

dada_denoise <- function(errorFile,derepFile,pairId,derep.rds=FALSE,paired=TRUE)
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

get_counts <- function(attribute) {
    get_counts_attr <- function(x) {
        #' Get the number of total and unique sequences
        #' on [attribute] attribute of object [x]
        total_count <- sum(x[[attribute]])
        uniq_count <- sum(x[[attribute]]>0)
        return(sprintf("%s (%s uniques)",total_count,uniq_count))
    }
    return(get_counts_attr)
}

write_merge_summary <- function(merged, key, samples) {
    ## Get a summary of matches
    data <- lapply(merged, function(df) table(df[[key]]))
    indices <- sort(unique(unlist(lapply(data, function(x) as.numeric(names(x))))))

    summary <- as.data.frame(matrix(0, nrow=length(data), ncol=1+length(indices)))
    summary[, 1] <- samples
    colnames(summary) <- c('sample', indices)

    for (i in 1:nrow(summary)) {
        summary[i, names(data[[i]])] <- as.numeric(data[[i]])
    }

    write.table(summary, sprintf("%s_summary.tsv", key), sep='\t', quote=FALSE, row.names=FALSE)
}

merge_reads <- function(minOverlap, maxMismatch)
{
    derepF <- list.files(path=".",pattern="*_R1.derep.RDS")
    derepR <- list.files(path=".",pattern="*_R2.derep.RDS")
    denoisedF <- list.files(path=".",pattern="*_R1.dada.RDS")
    denoisedR <- list.files(path=".",pattern="*_R2.dada.RDS")    
    
    inputs <- list(dadaF=lapply(denoisedF, readRDS),
                   derepF=lapply(derepF, readRDS),
                   dadaR=lapply(denoisedR, readRDS),
                   derepR=lapply(derepR, readRDS))
    params <- list(minOverlap=minOverlap,maxMismatch=maxMismatch,returnReject=TRUE)

    merged <- do.call(mergePairs,append(inputs,params))

    ## Get a summary of matches
    samples <- unname(sapply(denoisedF, function(x) gsub("_R1_dada.RDS", "", x)))
    write_merge_summary(merged, 'nmatch', samples)
    write_merge_summary(merged, 'nmismatch', samples)

    ## Save data
    saveRDS(merged,"reads_merged.RDS")

    ## Removes the reads not matching the merging criteria
    merged.accepted <- lapply(merged, function(df) df[df$accept,])

    return(merged.accepted)
}

make_summary <- function(steps, patterns, attributes) {
    sample.names <- as.character(sapply(
        list.files(pattern=patterns[1]),
        function(x) unlist(strsplit(x, gsub('\\*', '', patterns[1]),fixed=T))[1]
    ))
    summary <- data.frame(matrix(ncol=1+length(steps), nrow=length(sample.names)))
    colnames(summary) <- c("Sample",steps)
    summary['Sample'] = sample.names
    
    for( i in 1:length(patterns) ) {
        files <- list.files(pattern=patterns[i])
        data <- lapply(files, readRDS)
        summary[steps[i]] = sapply(data, get_counts(attributes[i]))
    }

    return(summary)
}
    
esv_table <- function(minOverlap=20, maxMismatch=1, paired=TRUE)
{
    if(paired) {
        merged <- merge_reads(minOverlap, maxMismatch)
    } else {
        ## Retrieve all the dada2 denoised files for merging
        denoisedFiles <- list.files(pattern="*_dada.RDS")
        pairIds <- as.character(sapply(denoisedFiles,
                                       function(x) unlist(strsplit(basename(x),"_"))[1]))
        merged <- lapply(denoisedFiles, function (x) readRDS(x))
        names(merged) <- pairIds
        saveRDS(merged, "fwd_reads_all.RDS")
    }

    esvTable <- makeSequenceTable(merged)

    ## Remove samples with no sequences after merging
    print("Removing empty samples")
    not_empty_samples <- rowSums(esvTable)>0
    esvTable <- esvTable[not_empty_samples, ]
    
    write.csv(esvTable,"raw_sequence_table.csv")
    
    esv_ids <- sprintf("esv_%s", c(1:dim(esvTable)[2]))
    uniquesToFasta(esvTable, "all_esv.fasta", ids=esv_ids)
    colnames(esvTable) <- esv_ids

    ## Make summary file
    if( length(list.files(pattern="*_nochimera.RDS")) > 0 ) {
        ## ITS pipeline
        steps = c("4-Dereplication","5-ChimeraRemoval","6-Denoising")
        patterns = c("*_R1.*_derep.RDS","*_R1.*_nochimera.RDS","*_R1.*_dada.RDS")
        attributes = c("uniques","uniques","denoised")
        summary <- make_summary(steps, patterns, attributes)

        write.csv(t(esvTable),"clustering_100.csv", quote=F)
    } else {
        ## 16S pipeline
        steps = c("3.1-Dereplication","3.2-Denoising")
        patterns = c("*_R1_derep.RDS","*_R1_dada.RDS")
        attributes = c("uniques","denoised")        
        summary <- make_summary(steps, patterns, attributes)

        if(paired) {
            summary['4-ESV'] = sapply(merged, get_counts('abundance'))
        } else {
            summary['4-ESV'] = summary['3.2-Denoising']
        }
        summary <- summary[not_empty_samples, ]
        
        ## write in .count_table format for Mothur
        esvTable <- cbind(esv_ids, colSums(esvTable), t(esvTable) )
        colnames(esvTable) <- c("Representative_Sequence","total",unlist(summary['Sample']))
        write.table(esvTable, file="all_esv.count_table", row.names=F, col.names=T, quote=F, sep="\t")
    }

    write.table(summary, "count_summary.tsv", row.names=F, sep="\t")                    
}

merge_fasta_for_vsearch <- function(filename)
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

lulu_curate <- function(abundanceFile,matchListFile,threshold,min_ratio_type,min_ratio,min_match,min_rel_cooccurence,cores=1)
{
    if(tools::file_ext(abundanceFile)=='shared'){
        otutab <- t(fast_table_load(abundanceFile, cores=cores,
                                    drop=c('label', 'numOtus'), row.label='Group'))
    } else {
        otutab <- fast_table_load(abundanceFile, cores=cores, row.label=1)        
    }
        
    abundance_table_name <- sprintf("all_lulu_%s.csv",threshold)
    curated_ids_name <- sprintf("lulu_ids_%s.csv",threshold)
    mapping_discarded_name <- sprintf("mapping_discarded_%s.txt", threshold)
    
    matchList <- read.table(matchListFile,
                            header=FALSE, col.names=c("OTU1","OTU2","pctIdentity"),
                            as.is=TRUE, check.names=F, stringsAsFactors=FALSE)

    if (dim(matchList)[1] > 0) {
        ## Run Lulu
        curated <- lulu(as.data.frame(otutab), matchList,
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

        write.table(curated$otu_map[curated$discarded_otus,],
                    mapping_discarded_name,
                    quote=F)

    } else {
        write.csv(otutab, abundance_table_name)

        write.table(rownames(otutab),
                    curated_ids_name,
                    quote=F,
                    row.names=F)

        empty.df <- read.table(text="", col.names=c('otu', 'parent_id'))
        write.table(empty.df, mapping_discarded_name)
    }

}

merge_otu_list <- function(list_file, merge_list, otu_thresh=100, cores=1) {
    list <- fast_table_load(list_file, row.label=-1, cores=cores)
    to_merge <- read.csv(merge_list, sep=' ', row.names=1, as.is=TRUE)

    for (daughter in rownames(to_merge)) {
        parent <- to_merge[daughter, 'parent_id']
        list[, parent] <- paste(list[, parent], list[, daughter], sep=',')
        list[, daughter] <- NULL
    }

    list[, 'numOtus'] <- dim(list)[2] - 2

    write.table(list, sprintf('all_lulu_%s.list', otu_thresh), sep='\t', quote=F, row.names=F)
}

get_species <- function(fasta_file, tax_file, db) {
    sequences <- read.fasta(fasta_file, 'DNA')
    sequences <- sapply(sequences, function(x) toupper(paste0(x[x!='-' & x!='.'], collapse='')))
    names(sequences) <- gsub('\t.*', '', names(sequences))

    tax <- read.table(tax_file, row.names=1, header=TRUE)$Size
    names(tax) <- sequences

    res <- assignSpecies(tax, db, allowMultiple=TRUE)
    res <- unname(res)

    rownames(res) <- names(sequences)
    colnames(res) <- c('Genus', 'Species')

    return(res)
}

calculate_unifrac <- function(tree_file, abundance_file, method='weighted', otu_thresh=100, cores=1) {

    registerDoParallel(cores=cores)

    abundance <- fast_table_load(abundance_file, cores=cores,
                                 drop=c('label', 'numOtus'), row.label='Group')
    tree <- read.tree(tree_file)

    ps <- phyloseq(
        otu_table(as.matrix(abundance), taxa_are_rows=FALSE),
        phy_tree(tree)
    )

    dists <- UniFrac(ps, weighted=(method=='weighted'), parallel=TRUE)

    write.csv(as.matrix(dists), sprintf('unifrac_%s_%s.csv', method, otu_thresh))
}
