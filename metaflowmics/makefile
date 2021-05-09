ifndef CONF
override CONF = local_test
endif

ifndef RELEASE
override RELEASE = latest
endif

test:
	nextflow run multithreaded-demux --inputdir "../tests/demux" -profile $(CONF) --n_per_file 100 --n_bases 1e3
	nextflow run Pipeline-ITS --reads "../tests/ITS/*R1.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF)

testITS:
	nextflow run Pipeline-ITS --reads "../tests/ITS/*R1.fastq.gz" -profile $(CONF) \
	&& nextflow run Pipeline-ITS --reads "../tests/ITS/*R{1,2}.fastq.gz" -profile $(CONF) --paired_end \
	&& nextflow run Pipeline-ITS --reads "../tests/ITS/*R1.fastq.gz" -profile $(CONF) --skip_lulu

test16S:
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF) --silva_db seed --compute_mothur_db \
	&& nextflow run Pipeline-16S --reads "../tests/16S/*_R1.fastq.gz" -profile $(CONF) --single_end --silva_db seed \
	&& nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF) --skip_subsampling --skip_lulu --silva_db seed \
	&& nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF) --skip_unifrac --silva_db seed

testdemux:
	nextflow run multithreaded-demux --inputdir "../tests/demux" -profile $(CONF) --n_per_file 100 --n_bases 1e3

clean:
	rm -rf work .nextflow*
	rm -rf demultiplexed
	rm -rf ITS-pipeline_outputs
	rm -rf 16S-pipeline_outputs

gcp_wrapper_16S:
	@ make gcp16S NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcp16S:
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp' --referenceAln databases/test.align --referenceTax databases/test.tax
gcp_wrapper_ITS:
	@ make gcpITS NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcpITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R1.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp'

setup_containers:
	module load tools/Singularity/3.7.1
	mkdir -p ~/.singularity_images.cache
	mkdir -p /tmp/sg && cd /tmp/sg
	singularity pull --docker-login docker://nakor/bash:5.1.4
	# docker://nakor/metaflowmics-r:0.0.1
	# docker://nakor/metaflowmics-python:0.0.1
	singularity pull quay.io/biocontainers/itsxpress:1.8.0--py_1
	# quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h399db7b_1
	# quay.io/biocontainers/mothur:1.44.1--hf0cea05_2
	# quay.io/biocontainers/vsearch:2.17.0--h95f258a_1
	mv *.img ~/.singularity_images.cache
