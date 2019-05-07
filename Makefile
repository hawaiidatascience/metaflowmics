test:
	nextflow run Pipeline-ITS -profile poire --reads "${PWD}/tests/ITS/*.fastq.gz"
	nextflow run Pipeline-16S -profile poire --reads "${PWD}/tests/16S/*_R{1,2}_*.fastq.gz"
testITS:
	nextflow run Pipeline-ITS -profile manoa_hpc --reads "${PWD}/tests/ITS/*.fastq.gz"
test16S:
	nextflow run Pipeline-16S -profile manoa_hpc --reads "${PWD}/tests/16S/*_R{1,2}_*.fastq.gz"
