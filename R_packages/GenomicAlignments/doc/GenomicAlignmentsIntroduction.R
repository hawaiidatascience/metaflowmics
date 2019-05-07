### R code from vignette source 'GenomicAlignmentsIntroduction.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## if (!require("BiocManager"))
##     install.packages("BiocManager")
## BiocManager::install("GenomicAlignments")


###################################################
### code chunk number 3: initialize
###################################################
library(GenomicAlignments)


###################################################
### code chunk number 4: readGAlignments
###################################################
library(GenomicAlignments)
aln1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
aln1 <- readGAlignments(aln1_file)
aln1
length(aln1)


###################################################
### code chunk number 5: accessors
###################################################
head(seqnames(aln1))
seqlevels(aln1)
head(strand(aln1))
head(cigar(aln1))
head(qwidth(aln1))
head(start(aln1))
head(end(aln1))
head(width(aln1))
head(njunc(aln1))


###################################################
### code chunk number 6: SessionInfo
###################################################
sessionInfo()


