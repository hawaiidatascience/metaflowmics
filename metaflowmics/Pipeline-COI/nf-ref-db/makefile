ROOT ?= $(PWD)/iBOL-database
ITER = 1
OUTDIR ?= $(ROOT)/iter-$(ITER)
DB ?= $(ROOT)/dl/cdhit/cdhit.repr.fa


dl:
	nextflow download_and_format.nf -w work/dl -profile local -resume \
	--outdir $(ROOT)/dl

align:
	nextflow progressive_alignment.nf -w work/updown -profile local -resume \
	--fasta $(FAA) \
	--outdir $(OUTDIR)

flag:
	python find_suspicious_aln.py \
	  --aln $(OUTDIR)/down/species.updated.afa \
	  --outdir $(OUTDIR)

dna_align:
	sed '/^>/!s/[-.]//g' $(OUTDIR)/down/species.updated.afa > $(OUTDIR)/species.faa
	cd ../../subworkflows \
	&& nextflow hmmer.nf -w $(PWD)/work/dna_aln -profile local -resume \
	--fna $(ROOT)/dl/iBOL_COI.fna \
	--faa $(OUTDIR)/species.faa \
	--db $(OUTDIR)/down/species.updated.afa \
	--outdir $(OUTDIR)/MSA_backtranslated \
	&& rm -rf ../../subworkflows/.nextflow*

all:
	@make dl
	@make align ITER=1 FAA=$(ROOT)/dl/cdhit/cdhit.repr.fa
	@make flag ITER=1
	@make align ITER=2 FAA=$(ROOT)/iter-1/species.updated.kept.faa
	@make dna_align ITER=2

mothur:
	sed 's/,/;/g' $(AFA) | sed -r 's/[kpcofgs]f?__//g' | sed -r 's/>(.*)(id__[^ ]*)(.*)/>\2\t\3\t\1/g' > $(AFA).align
	grep '^>' $(AFA).align | cut -c2- | cut -f1,3 > $(AFA).tax

