test:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF)
testITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R{1,2}.fastq.gz" --pairedEnd -profile $(CONF)
test16S:
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" --skipSubsampling -profile $(CONF)

python_container:
	sudo docker build -f Dockerfile_python -t nakor/python_libs_pipeline scripts/
	sudo docker push nakor/python_libs_pipeline
R_container:
	sudo docker build -f Dockerfile_R -t nakor/r_libs_pipeline scripts/
	sudo docker push nakor/r_libs_pipeline
mothur_container:
	sudo docker build -f Dockerfile_mothur -t nakor/mothur_pipeline scripts/
	sudo docker push nakor/mothur_pipeline
gcp_wrapper_16S:
	@ make gcp16S NXF_VER=19.01.0 NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcp16S:
	nextflow1901 run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp'
gcp_wrapper_ITS:
	@ make gcpITS NXF_VER=19.01.0 NXF_MODE=google GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp_credentials.json"
gcpITS:
	nextflow1901 run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R1.fastq.gz" -profile gcp -work-dir 'gs://c-maiki-work/tmp'
