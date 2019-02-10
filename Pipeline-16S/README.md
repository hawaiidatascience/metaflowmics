# 16S-pipeline

## Getting started

Follow these instructions to get the pipeline started on your machine

### Pre-requisites

To run the pipeline, you will need to satisfy the following dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- python3 + libraries: Biopython, pandas
- R + libraries: ggplot2, lulu, dada2, seqinr, stringr
- [Mothur](https://github.com/mothur/mothur) (tested with v1.41.2) 

### Usage

- Clone the repository:
`git clone https://github.com/Puumanamana/16S-pipeline.git`

To run the pipeline locally, you need to set up a configuration file. An example is available in `conf/poire.config`
Pipeline parameters can be set in the file `nextflow.config`

Then, you can run the pipeline by running:
`nextflow run 16S-pipeline -profile manoa_hpc --revRead 1 --reads PATH_TO_READS/PATTERN`

## Pipeline summary

**dada2**
- filterAndTrim (minimum read length, maximum expected errors, read length truncation, quality filter)
- learnErrors
- mergePairs (minimum overlap and maximum mismatches)

**Mothur**
- align.seqs; filter.seqs; screen.seqs  --> alignment against reference for further filtering
- chimera.vsearch; remove.seqs          --> chimera removal
- sub.sample
- cluster (at 95,97,99 and 100% identity)
- classify.seqs ; classify.otu 

**Lulu**
- preLulu (create matchlists for LULU)
- LULU

**Python**
- OutputFilter (remove LULU flagged sequences from FASTA and taxonomy file; remove low abundant OTUs (<=1))

**Postprocessing**
- ConvertToMothur: Convert output file to .shared file
- Results: Mothur postprocessing with get.relabund, clearcut and unifrac.weighted
- SummaryFile: Generates summary of reads per sample per step
