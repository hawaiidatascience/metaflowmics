ROOT ?= $(PWD)/iBOL-database
ITER = 2
OUTDIR ?= $(ROOT)/iter-$(ITER)
DB ?= $(ROOT)/dl/cdhit/cdhit.repr.fa


dl:
	nextflow download_and_format.nf -w work/dl -profile local -resume \
	--outdir $(ROOT)/dl $(ARGS)

rdp:
	grep '>' $(ROOT)/iBOL_COI.fna | sed 's/ /;/' | cut -d' ' -f1 | tr ';' $'\t' | cut -c2- > rdp/taxonomy.tsv
	python py/fix_tax.py --tax rdp/taxonomy.tsv
	sed 's/immidjanzen01/immidJanzen01/g' rdp/taxonomy_clean.tsv | sed 's/Enchytraeidaegen/EnchytraeidaeGEN/g' | sed 's/oecophor01/Oecophor01/g' > tmp && mv -f tmp rdp/taxonomy_clean.tsv
	python filter_fasta.py --fasta $(ROOT)/iBOL_COI.fna --tax rdp/taxonomy_clean.tsv
	python rdp_py_scripts/addFullLineage.py rdp/taxonomy_clean.tsv rdp/rawSeq.fasta > rdp/ready4train_seqs.fasta
	python rdp_py_scripts/lineage2taxTrain.py rdp/taxonomy_clean.tsv > rdp/ready4train_taxonomy.txt

	rdp_classifier -Xmx256g train -o rdp/rdp_COI_trained -s rdp/ready4train_seqs.fasta -t rdp/ready4train_taxonomy.txt

align:
	nextflow progressive_alignment.nf -w work/updown -profile local -resume \
	--fasta $(FAA) \
	--outdir $(OUTDIR)

flag:
	python py/find_suspicious_aln.py \
	  --aln $(OUTDIR)/down/species.updated.afa \
	  --outdir $(OUTDIR)

dna_align:
	sed '/^>/!s/[-.]//g' $(OUTDIR)/down/species.updated.afa > $(OUTDIR)/species.faa \
	&& python py/synchronize_fasta.py \
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
	sed -r 's/[kpcofgs]f?__//g' $(INDIR)/species.codons.afa | \
	    sed -r 's/>(.*);(id__[^ ]*)(.*)/>\2\t100\t\1;/g' > $(INDIR)/db.codons.align \
	&& grep '^>' $(INDIR)/db.codons.align | cut -c2- | cut -f1,3 \
	> $(INDIR)/db.codons.tax

all:
	@make dl ARGS="--taxa $(HOME)/db/BOLD/taxa_to_dl.txt --external $(PWD)/manual"
	@make align ITER=1 FAA=$(ROOT)/dl/cdhit/cdhit.repr.fa
	@make flag ITER=1
	@make align ITER=2 FAA=$(ROOT)/iter-1/species.updated.kept.faa
	@make dna_align ITER=2
	@make mothur INDIR=$(ROOT)/iter-2/MSA_backtranslated/tranalign

