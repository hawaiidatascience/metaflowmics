#!/usr/bin/env bash

test:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile manoa_hpc 
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile manoa_hpc
testITS:
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*R1.fastq.gz" -profile manoa_hpc
	nextflow run Pipeline-ITS --reads "${PWD}/tests/ITS/*_R{1,2}.fastq.gz" --pairedEnd -profile manoa_hpc
test16S:
	nextflow run Pipeline-16S --reads "${PWD}/tests/16S/*_R{1,2}.fastq.gz" -profile manoa_hpc
