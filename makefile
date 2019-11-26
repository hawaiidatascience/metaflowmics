test:
	nextflow run multithreaded_demux --inputdir "${PWD}/tests/demux" -profile $(CONF) --n_per_file 100 --n_bases 1e3
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" --referenceAln databases/silva.seed_v132.align --referenceTax databases/silva.seed_v132.tax -profile $(CONF)

testITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R{1,2}.fastq.gz" --pairedEnd -profile $(CONF)

test16S:
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" --referenceAln databases/silva.seed_v132.align --referenceTax databases/silva.seed_v132.tax -profile $(CONF)
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R1.fastq.gz" --referenceAln databases/silva.seed_v132.align --referenceTax databases/silva.seed_v132.tax --singleEnd --truncLen 240 -profile $(CONF)

testdemux:
	nextflow run multithreaded_demux --inputdir "${PWD}/tests/demux" -profile $(CONF) --n_per_file 100 --n_bases 1e3

clean:
	rm -rf work .nextflow*
	rm -rf demultiplexed
	rm -rf ITS-pipeline_outputs
	rm -rf ITS-pipeline_outputs

python_container:
	sudo docker build -f python.dockerfile -t nakor/python_libs_pipeline scripts/
	sudo docker push nakor/python_libs_pipeline
R_container:
	sudo docker build -f r.dockerfile -t nakor/r_libs_pipeline scripts/
	sudo docker push nakor/r_libs_pipeline
mothur_container:
	sudo docker build -f mothur.dockerfile -t nakor/mothur_pipeline scripts/
	sudo docker push nakor/mothur_pipeline

gcp_wrapper_16S:
	@ make gcp16S NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcp16S:
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp' --referenceAln databases/silva.seed_v132.align --referenceTax databases/silva.seed_v132.tax
gcp_wrapper_ITS:
	@ make gcpITS NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcpITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R1.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp'
