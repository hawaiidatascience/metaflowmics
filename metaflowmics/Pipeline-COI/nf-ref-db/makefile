OUTDIR ?= db_COI
DB ?= $(OUTDIR)/dl/cdhit/test.repr.fa

dl:
	nextflow download_and_format.nf -w work/dl -profile local -resume \
	--outdir $(OUTDIR)/dl

main:
	nextflow main.nf -w work/updown -profile local -resume \
	--fasta $(DB) \
	--outdir $(OUTDIR)
