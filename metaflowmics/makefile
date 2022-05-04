CONF ?= local_test
PATTERN ?= R{1,2}
ARGS ?=

test:
	@make testITS PATTERN=R1
	@make test16S
	@make testdemux
	@make testCOI

test16S:
	mkdir -p testrun
	cd $(shell mktemp -d -p testrun) && \
	nextflow run $(PWD)/Pipeline-16S -resume -profile $(CONF) \
	  --reads "$(PWD)/../tests/16S/*_$(PATTERN).fastq.gz" --silva_release seed $(ARGS)

testITS:
	mkdir -p testrun
	cd $(shell mktemp -d -p testrun) && \
	 nextflow run $(PWD)/Pipeline-ITS -resume -profile $(CONF) \
	  --reads "$(PWD)/../tests/ITS/*_$(PATTERN).fastq.gz" $(ARGS)

testCOI:
	mkdir -p testrun
	cd $(shell mktemp -d -p testrun) && \
	nextflow run $(PWD)/Pipeline-COI -resume -profile $(CONF) \
	  --reads "$(PWD)/../tests/COI/*_$(PATTERN).fastq.gz" \
	  --barcodes $(PWD)/../tests/COI/index.fasta \
	  --db_dir $(PWD)/../db/COI \
	  --min_subsampling 100 --demux jamp $(ARGS)

testdemux:
	mkdir -p testrun
	cd $(shell mktemp -d -p testrun) && \
	nextflow run $(PWD)/multithreaded-demux -resume -profile $(CONF) \
	  --inputdir "$(PWD)/../tests/demux" --n_per_file 100 --n_bases 1e3 $(ARGS)

test16Sall:
	@make test16S
	@make test16S ARGS="--pool T --compute_mothur_db --skip_unifrac"
	@make test16S ARGS="--single_end --skip_subsampling --skip_lulu"

testITSall:
	@make testITS
	@make testITS ARGS="--pool T --paired_end --skip_lulu"
	@make testITS ARGS="--skip_lulu"

clean:
	rm -rf work .nextflow* testrun

gcp_wrapper_16S:
	@ make gcp16S NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"

gcp16S:
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp' --referenceAln databases/test.align --referenceTax databases/test.tax

gcp_wrapper_ITS:
	@ make gcpITS NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcpITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R1.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp'
