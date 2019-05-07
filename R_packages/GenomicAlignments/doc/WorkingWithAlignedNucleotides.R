### R code from vignette source 'WorkingWithAlignedNucleotides.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: bamfiles
###################################################
library(RNAseqData.HNRNPC.bam.chr14)
bamfiles <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
names(bamfiles)  # the names of the runs


###################################################
### code chunk number 3: quickBamFlagSummary
###################################################
library(Rsamtools)
quickBamFlagSummary(bamfiles[1], main.groups.only=TRUE)


###################################################
### code chunk number 4: ScanBamParam
###################################################
flag1 <- scanBamFlag(isFirstMateRead=TRUE, isSecondMateRead=FALSE,
                     isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
param1 <- ScanBamParam(flag=flag1, what="seq")


###################################################
### code chunk number 5: readGAlignments
###################################################
library(GenomicAlignments)
gal1 <- readGAlignments(bamfiles[1], use.names=TRUE, param=param1)


###################################################
### code chunk number 6: read_sequences
###################################################
mcols(gal1)$seq


###################################################
### code chunk number 7: original-query-sequences
###################################################
oqseq1 <- mcols(gal1)$seq
is_on_minus <- as.logical(strand(gal1) == "-")
oqseq1[is_on_minus] <- reverseComplement(oqseq1[is_on_minus])


###################################################
### code chunk number 8: is_dup
###################################################
is_dup <- duplicated(names(gal1))
table(is_dup)


###################################################
### code chunk number 9: same-name-implies-same-seq-in-U1-oqseq
###################################################
dup2unq <- match(names(gal1), names(gal1))
stopifnot(all(oqseq1 == oqseq1[dup2unq]))


###################################################
### code chunk number 10: oqseq1
###################################################
oqseq1 <- oqseq1[!is_dup]


###################################################
### code chunk number 11: most_frequent_cigars
###################################################
head(sort(table(cigar(gal1)), decreasing=TRUE))


###################################################
### code chunk number 12: cigarOpTable
###################################################
colSums(cigarOpTable(cigar(gal1)))


###################################################
### code chunk number 13: table_njunc
###################################################
table(njunc(gal1))


###################################################
### code chunk number 14: readGAlignmentPairs
###################################################
library(pasillaBamSubset)
flag0 <- scanBamFlag(isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
param0 <- ScanBamParam(flag=flag0)
U3.galp <- readGAlignmentPairs(untreated3_chr4(), use.names=TRUE, param=param0)
head(U3.galp)


###################################################
### code chunk number 15: first-and-last-U3.galp
###################################################
head(first(U3.galp))
head(last(U3.galp))


###################################################
### code chunk number 16: isProperPair
###################################################
table(isProperPair(U3.galp))


###################################################
### code chunk number 17: keep-only-proper-pairs
###################################################
U3.GALP <- U3.galp[isProperPair(U3.galp)]


###################################################
### code chunk number 18: U3.GALP_names_is_dup
###################################################
U3.GALP_names_is_dup <- duplicated(names(U3.GALP))
table(U3.GALP_names_is_dup)


###################################################
### code chunk number 19: U3.GALP_qnames
###################################################
U3.uqnames <- unique(names(U3.GALP))
U3.GALP_qnames <- factor(names(U3.GALP), levels=U3.uqnames)


###################################################
### code chunk number 20: U3.GALP_dup2unq
###################################################
U3.GALP_dup2unq <- match(U3.GALP_qnames, U3.GALP_qnames)


###################################################
### code chunk number 21: gaps-in-U3.GALP
###################################################
head(unique(cigar(first(U3.GALP))))
head(unique(cigar(last(U3.GALP))))
table(njunc(first(U3.GALP)), njunc(last(U3.GALP)))


###################################################
### code chunk number 22: no-indels-in-U3.GALP
###################################################
colSums(cigarOpTable(cigar(first(U3.GALP))))
colSums(cigarOpTable(cigar(last(U3.GALP))))


###################################################
### code chunk number 23: sessionInfo
###################################################
sessionInfo()


