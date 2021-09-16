ROOT ?= $(PWD)/iBOL-database
ITER = 2
OUTDIR ?= $(ROOT)/iter-$(ITER)
DB ?= $(ROOT)/dl/cdhit/cdhit.repr.fa


dl:
	nextflow download_and_format.nf -w work/dl -profile local -resume \
	--outdir $(ROOT)/dl $(ARGS)

align:
	nextflow progressive_alignment.nf -w work/updown -profile local -resume \
	--fasta $(FAA) \
	--outdir $(OUTDIR)

flag:
	python find_suspicious_aln.py \
	  --aln $(OUTDIR)/down/species.updated.afa \
	  --outdir $(OUTDIR)

dna_align:
	sed '/^>/!s/[-.]//g' $(OUTDIR)/down/species.updated.afa > $(OUTDIR)/species.faa \
	&& python synchronize_fasta.py \
	--query $(ROOT)/dl/iBOL_COI.fna \
	--ref $(OUTDIR)/species.faa \
	--output $(OUTDIR)/species.fna \
	&& cd ../../subworkflows \
	&& nextflow hmmer.nf -w $(PWD)/work/dna_aln -profile local -resume \
	--fna $(OUTDIR)/species.fna \
	--faa $(OUTDIR)/species.faa \
	--db $(OUTDIR)/down/species.updated.afa \
	--outdir $(OUTDIR)/MSA_backtranslated \
	&& rm -rf ../../subworkflows/.nextflow*

mothur:
	./convert_to_mothur.sh $(INDIR)/species.codons.afa

all:
	@make dl
	@make align ITER=1 FAA=$(ROOT)/dl/cdhit/cdhit.repr.fa
	@make flag ITER=1
	@make align ITER=2 FAA=$(ROOT)/iter-1/species.updated.kept.faa
	@make dna_align ITER=2
	@make mothur INDIR=$(ROOT)/iter-2/MSA_backtranslated/tranalign

