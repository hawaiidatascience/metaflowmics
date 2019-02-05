# ITS-pipeline

## Getting started

Follow these instructions to get the pipeline started on your machine

### Pre-requisites

To run the pipeline, you will need to satisfy the following dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- python3 + libraries: Biopython, pandas
- R + libraries: ggplot2, lulu, dada2, seqinr, stringr
- Mothur (tested with v1.41.2) 

### Usage

- Clone the repository:
`git clone https://github.com/Puumanamana/16S-pipeline.git`

To run the pipeline locally, you need to set up a configuration file. An example is available in `conf/poire.config`
Pipeline parameters can be set in the file `nextflow.config`

Then, you can run the pipeline by running:
`nextflow run 16S-pipeline -profile manoa_hpc --revRead 1 --reads PATH_TO_READS/PATTERN`

