CONF ?= local,singularity

test:
	mkdir -p testrun
	nextflow run Pipeline-16S -w testrun/work -resume -profile $(CONF),test16S --outdir testrun/16S
	nextflow run Pipeline-ITS -w testrun/work -resume -profile $(CONF),testITS --outdir testrun/ITS
	nextflow run Pipeline-COI -w testrun/work -resume -profile $(CONF),testCOI --outdir testrun/COI
	nextflow run multithreaded-demux -w testrun/work -resume -profile $(CONF),testDemux --outdir testrun/demux

clean:
	rm -rf .nextflow* testrun

gcp_wrapper_16S:
	@ make gcp16S NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"

gcp16S:
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp' --referenceAln databases/test.align --referenceTax databases/test.tax

gcp_wrapper_ITS:
	@ make gcpITS NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcpITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R1.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp'
