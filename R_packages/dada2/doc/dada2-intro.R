## ----filenames, message=FALSE, warning=FALSE-------------------------------
library(dada2); packageVersion("dada2")
fnF1 <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
fnR1 <- system.file("extdata", "sam1R.fastq.gz", package="dada2")
filtF1 <- tempfile(fileext=".fastq.gz")
filtR1 <- tempfile(fileext=".fastq.gz")

## ----inspect---------------------------------------------------------------
plotQualityProfile(fnF1) # Forward
plotQualityProfile(fnR1) # Reverse

## ----filter----------------------------------------------------------------
filterAndTrim(fwd=fnF1, filt=filtF1, rev=fnR1, filt.rev=filtR1,
                  trimLeft=10, truncLen=c(240, 200), 
                  maxN=0, maxEE=2,
                  compress=TRUE, verbose=TRUE)

## ----derep-----------------------------------------------------------------
derepF1 <- derepFastq(filtF1, verbose=TRUE)
derepR1 <- derepFastq(filtR1, verbose=TRUE)

## ----learn, warning=FALSE--------------------------------------------------
errF <- learnErrors(derepF1, multithread=FALSE) # multithreading is available on many functions
errR <- learnErrors(derepR1, multithread=FALSE)

## ----dada, warning=FALSE---------------------------------------------------
dadaF1 <- dada(derepF1, err=errF, multithread=FALSE)
dadaR1 <- dada(derepR1, err=errR, multithread=FALSE)
print(dadaF1)

## ----merge-----------------------------------------------------------------
merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose=TRUE)

## ----bimeras---------------------------------------------------------------
merger1.nochim <- removeBimeraDenovo(merger1, multithread=FALSE, verbose=TRUE)

## ----sample2, warning=FALSE------------------------------------------------
# Assign filenames
fnF2 <- system.file("extdata", "sam2F.fastq.gz", package="dada2")
fnR2 <- system.file("extdata", "sam2R.fastq.gz", package="dada2")
filtF2 <- tempfile(fileext=".fastq.gz")
filtR2 <- tempfile(fileext=".fastq.gz")
# Filter and Trim
filterAndTrim(fwd=fnF2, filt=filtF2, rev=fnR2, filt.rev=filtR2, maxN=0, trimLeft=10, truncLen=c(240, 200), maxEE=2, compress=TRUE, verbose=TRUE)
# Dereplicate
derepF2 <- derepFastq(filtF2, verbose=TRUE)
derepR2 <- derepFastq(filtR2, verbose=TRUE)
# Infer sample composition (using already learned error rates)
dadaF2 <- dada(derepF2, err=errF, multithread=FALSE)
dadaR2 <- dada(derepR2, err=errR, multithread=FALSE)
# Merge
merger2 <- mergePairs(dadaF2, derepF2, dadaR2, derepR2, verbose=TRUE)

## ----make-table------------------------------------------------------------
seqtab <- makeSequenceTable(list(merger1, merger2))
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)

