library(ShortRead)

denoised2fastq <- function(filename, output) {
    dada.obj <- readRDS(filename)

    ids <- paste0('seq', 1:length(dada.obj$denoised))
    sequences <- dada.obj$sequence
    qualities <- apply(dada.obj$quality, 1, function(x) intToUtf8(x+33))

    reads <- ShortReadQ(
        DNAStringSet(sequences),
        FastqQuality(qualities),
        BStringSet(ids)
    )

    writeFastq(reads, output, compress=F)
}

dada2fastq <- function(derepF.file, derepR.file, dadaF.file, dadaR.file, sample) {
    derep <- list(fwd=readRDS(derepF.file), rev=readRDS(derepR.file))
    denoised <- list(fwd=readRDS(dadaF.file), rev=readRDS(dadaR.file))

    pairdf <- data.frame(fwd=denoised$fwd$map[derep$fwd$map],
                         rev=denoised$rev$map[derep$rev$map])

    pairdf$rids <- sprintf("%s__%s", sample, 1:dim(pairdf)[1])

    pairdf <- pairdf[!is.na(pairdf$fwd) & !is.na(pairdf$rev), ]

    for(key in c('fwd', 'rev')) {
        nb <- as.integer(key=='fwd') + 1

        repr_ids <- pairdf[[key]]
        
        ids <- sprintf('%s %s:N:0:0', pairdf$rid, nb)
        sequences <- denoised[[key]]$sequence[repr_ids]
        qualities <- apply(denoised[[key]]$quality[repr_ids, ], 1, function(x) intToUtf8(x+33))

        reads <- ShortReadQ(
            DNAStringSet(sequences),
            FastqQuality(qualities),
            BStringSet(ids)
        )
        
        writeFastq(reads, sprintf('%s_R%s_denoised.fastq', sample, nb), compress=F)
    }
}

convert_all <- function() {
    sample.names <- as.character(sapply(
        list.files(path=".",pattern="*_R1_derep.RDS"),
        function(x) unlist(strsplit(x, "_R1_derep.RDS", fixed=T))[1]))

    for(sample in sample.names) {
        derepF <- sprintf("%s_R1_derep.RDS", sample)
        derepR <- sprintf("%s_R2_derep.RDS", sample)
        dadaF <- sprintf("%s_R1_dada.RDS", sample)
        dadaR <- sprintf("%s_R2_dada.RDS", sample)                
        dada2fastq(derepF, derepR, dadaF, dadaR, sample)
    }

    info.file <- data.frame(
        sample=sample.names,
        fwd=sprintf('%s_R1_denoised.fastq', sample.names),
        rev=sprintf('%s_R2_denoised.fastq', sample.names)
    )

    write.table(info.file, 'info.file', col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
}
