# ITS-pipeline

## Getting started

Follow these instructions to get the pipeline started on your machine

### Pre-requisites

To run the pipeline, you will need to satisfy the following dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- python3 + libraries: Biopython, pandas
- R + libraries: ggplot2, lulu, dada2, seqinr, stringr
- [FastQC](https://github.com/s-andrews/FastQC)
- [ITSxpress](https://github.com/USDA-ARS-GBRU/itsxpress) (tested with v1.7.2)
- [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/download.html) (tested with v0.0.14)
- [vsearch](https://github.com/torognes/vsearch) (tested with v2.10.4)
- The modified CONSTAX algorithm

The CONSTAX algorithm requires a few configuration steps:
- Set up usearch8 and usearch10 available [here](https://www.drive5.com/usearch/download.html)    
- RDPTools (directions available [here](https://github.com/natalie-vandepol/compare_taxonomy/tree/master/CONSTAX))

Once installed, modify the configuration file in friendly_CONSTAX/constax.sh with:
- RDPTools_path: path to the RDPTools folder

### Usage
`nextflow run ITS-pipeline -profile manoa_hpc --its ITS1 --revRead 0` (or 1 to include the reverse reads)

