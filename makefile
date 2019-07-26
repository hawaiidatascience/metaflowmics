test:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF)
testITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile $(CONF)
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R{1,2}.fastq.gz" --pairedEnd -profile $(CONF)
test16S:
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile $(CONF)

python_container:
	sudo docker build -f Dockerfile_python -t nakor/python_libs_pipeline Pipeline-16S/scripts/
	sudo docker push nakor/python_libs_pipeline
R_container:
	sudo docker build -f Dockerfile_R -t nakor/r_libs_pipeline Pipeline-16S/scripts/
	sudo docker push nakor/r_libs_pipeline
mothur_container:
	sudo docker build -f Dockerfile_mothur -t nakor/mothur_pipeline Pipeline-16S/scripts/
	sudo docker push nakor/mothur_pipeline
