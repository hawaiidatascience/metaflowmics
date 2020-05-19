ifndef CONF
override CONF = local_test
endif

ifndef RELEASE
override RELEASE = latest
endif

test:
	nextflow run multithreaded-demux --inputdir "../tests/demux" -profile $(CONF) --n_per_file 100 --n_bases 1e3
	nextflow run Pipeline-ITS --reads "../tests/ITS/*R1.fastq.gz" -profile $(CONF) --uniteDB databases/test/test_unite.fasta
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" --referenceAln databases/test/test.align --referenceTax databases/test/test.tax -profile $(CONF)

testITS:
	nextflow run Pipeline-ITS --reads "../tests/ITS/*R1.fastq.gz" -profile $(CONF) --uniteDB databases/test/test_unite.fasta
	# nextflow run Pipeline-ITS --reads "../tests/ITS/*_R{1,2}.fastq.gz" --pairedEnd -profile $(CONF) --uniteDB databases/test/test_unite.fasta

test16S:
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" --referenceAln databases/test/test.align --referenceTax databases/test/test.tax -profile $(CONF)
	# nextflow run Pipeline-16S --reads "../tests/16S/*_R1.fastq.gz" --referenceAln databases/test/test.align --referenceTax databases/test/test.tax --singleEnd --truncLen 240 -profile $(CONF)

testdemux:
	nextflow run multithreaded-demux --inputdir "../tests/demux" -profile $(CONF) --n_per_file 100 --n_bases 1e3

clean:
	rm -rf work .nextflow*
	rm -rf demultiplexed
	rm -rf ITS-pipeline_outputs
	rm -rf 16S-pipeline_outputs

container:
	sudo docker build -f ../pipeline-env.dockerfile -t nakor/pipeline-env:$(RELEASE) scripts/
	sudo docker push nakor/pipeline-env:$(RELEASE)
	rm -f ~/.singularity_images.cache/nakor-pipeline-env.img

gcp_wrapper_16S:
	@ make gcp16S NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcp16S:
	nextflow run Pipeline-16S --reads "../tests/16S/*_R{1,2}.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp' --referenceAln databases/test.align --referenceTax databases/test.tax
gcp_wrapper_ITS:
	@ make gcpITS NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcpITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R1.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp' --uniteDB databases/test_unite.fasta
