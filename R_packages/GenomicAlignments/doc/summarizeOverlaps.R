### R code from vignette source 'summarizeOverlaps.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: options
###################################################
options(width=72)
options("showHeadLines" = 3)
options("showTailLines" = 3)


###################################################
### code chunk number 3: firstExample
###################################################
library(GenomicAlignments)
library(DESeq2)
library(edgeR)

fls <- list.files(system.file("extdata", package="GenomicAlignments"),
    recursive=TRUE, pattern="*bam$", full=TRUE)

features <- GRanges(
    seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
    ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 4000, 
        7500, 5000, 5400), width=c(rep(500, 3), 600, 900, 500, 300, 900, 
        300, 500, 500)), "-",
    group_id=c(rep("A", 4), rep("B", 5), rep("C", 2)))

olap <- summarizeOverlaps(features, fls)
deseq <- DESeqDataSet(olap, design= ~ 1)
edger <- DGEList(assay(olap), group=rownames(colData(olap)))


###################################################
### code chunk number 4: simple
###################################################
rd <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(100),
    cigar = "300M", strand = strand("+"))

gr1 <- GRanges("chr1", IRanges(start=50, width=150), strand="+")
gr2 <- GRanges("chr1", IRanges(start=350, width=150), strand="+")


###################################################
### code chunk number 5: simpleGRanges
###################################################
gr <- c(gr1, gr2)
data.frame(union = assay(summarizeOverlaps(gr, rd)),
           intStrict = assay(summarizeOverlaps(gr, rd,
               mode="IntersectionStrict")),
           intNotEmpty = assay(summarizeOverlaps(gr, rd,
               mode="IntersectionNotEmpty")))


###################################################
### code chunk number 6: simpleGRangesList
###################################################
grl <- GRangesList(c(gr1, gr2))
data.frame(union = assay(summarizeOverlaps(grl, rd)),
           intStrict = assay(summarizeOverlaps(grl, rd,
               mode="IntersectionStrict")),
           intNotEmpty = assay(summarizeOverlaps(grl, rd,
               mode="IntersectionNotEmpty")))


###################################################
### code chunk number 7: data
###################################################
group_id <- c("A", "B", "C", "C", "D", "D", "E", "F", "G", "G", "H", "H")
features <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr1", "chr2", "chr2",
        "chr1", "chr1", "chr2", "chr2", "chr1", "chr1")),
    strand = strand(rep("+", length(group_id))),
    ranges = IRanges(
        start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000, 3000, 3350, 5000, 5400),
        width=c(500, 900, 500, 300, 600, 300, 500, 900, 150, 200, 500, 500)),
   DataFrame(group_id)
)

reads <- GAlignments(
    names = c("a","b","c","d","e","f","g"),
    seqnames = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
    pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    cigar = c("500M", "100M", "300M", "500M", "300M", "50M200N50M", "50M150N50M"),
    strand = strand(rep.int("+", 7L)))



###################################################
### code chunk number 8: GRanges
###################################################
data.frame(union = assay(summarizeOverlaps(features, reads)),
           intStrict = assay(summarizeOverlaps(features, reads,
               mode="IntersectionStrict")),
           intNotEmpty = assay(summarizeOverlaps(features, reads,
               mode="IntersectionNotEmpty")))


###################################################
### code chunk number 9: lst
###################################################
lst <- split(features, mcols(features)[["group_id"]])
length(lst)


###################################################
### code chunk number 10: GRangesList
###################################################
data.frame(union = assay(summarizeOverlaps(lst, reads)),
           intStrict = assay(summarizeOverlaps(lst, reads,
               mode="IntersectionStrict")),
           intNotEmpty = assay(summarizeOverlaps(lst, reads,
               mode="IntersectionNotEmpty")))


###################################################
### code chunk number 11: gff (eval = FALSE)
###################################################
## library(rtracklayer)
## fl <- paste0("ftp://ftp.ensembl.org/pub/release-62/",
##              "gtf/drosophila_melanogaster/",
##              "Drosophila_melanogaster.BDGP5.25.62.gtf.gz")
## gffFile <- file.path(tempdir(), basename(fl))
## download.file(fl, gffFile)
## gff0 <- import(gffFile)


###################################################
### code chunk number 12: gff_parse (eval = FALSE)
###################################################
## idx <- mcols(gff0)$source == "protein_coding" & 
##            mcols(gff0)$type == "exon" & 
##            seqnames(gff0) == "4"
## gff <- gff0[idx]
## ## adjust seqnames to match Bam files
## seqlevels(gff) <- paste("chr", seqlevels(gff), sep="")
## chr4genes <- split(gff, mcols(gff)$gene_id)


###################################################
### code chunk number 13: pasilla_param
###################################################
param <- ScanBamParam(
             what='qual',
             which=GRanges("chr4", IRanges(1, 1e6)),
             flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=NA),
             tag="NH")


###################################################
### code chunk number 14: pasilla_count (eval = FALSE)
###################################################
## fls <- c("treated1.bam", "untreated1.bam", "untreated2.bam")
## path <- "pathToBAMFiles"
## bamlst <- BamFileList(fls)
## genehits <- summarizeOverlaps(chr4genes, bamlst, mode="Union")


###################################################
### code chunk number 15: pasilla_exoncountset (eval = FALSE)
###################################################
## expdata <- MIAME(
##               name="pasilla knockdown",
##               lab="Genetics and Developmental Biology, University of 
##                   Connecticut Health Center",
##               contact="Dr. Brenton Graveley",
##               title="modENCODE Drosophila pasilla RNA Binding Protein RNAi 
##                   knockdown RNA-Seq Studies",
##               pubMedIds="20921232",
##               url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE18508",
##               abstract="RNA-seq of 3 biological replicates of from the Drosophila
##                   melanogaster S2-DRSC cells that have been RNAi depleted of mRNAs 
##                   encoding pasilla, a mRNA binding protein and 4 biological replicates 
##                   of the the untreated cell line.")
## 
## design <- data.frame(
##               condition=c("treated", "untreated", "untreated"),
##               replicate=c(1,1,2),
##               type=rep("single-read", 3),
##               countfiles=path(colData(genehits)[,1]), stringsAsFactors=TRUE)
## 
## geneCDS <- DESeqDataSet(genehits, design=design, metadata=list(expdata=expdata))


###################################################
### code chunk number 16: pasilla_genes (eval = FALSE)
###################################################
## chr4tx <- split(gff, mcols(gff)$transcript_id)
## txhits <- summarizeOverlaps(chr4tx, bamlst)
## txCDS <- DESeqDataSet(txhits, design=design, metadata=list(expdata=expdata))


