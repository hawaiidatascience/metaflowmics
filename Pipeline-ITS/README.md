# ITS-pipeline

## Getting started

Follow these instructions to get the pipeline started on your machine

### Pre-requisites

To run the pipeline, you will need to satisfy the following dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- python3 + libraries: Biopython, pandas
- R + libraries: ggplot2, lulu, dada2, seqinr, stringr, ShortRead
- [Mothur](https://github.com/mothur/mothur) (tested with v1.43) 

### Usage

- Clone the repository:
```
git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
cd Pipeline-ITS
```

#### Make your own configuration file
To run the pipeline locally, you need to set up a configuration file. An example is available in `conf/local.config`.
You can modify this configuration to fit the specs of your machine.

#### Changing the default pipeline parameters

Pipeline parameters can be set either:
- in the file `nextflow.config`
- by using the flags --[PARAMETER-NAME] [PARAMETER-VALUE] (see `nextflow run ITS-pipeline --help`)

### Running the pipeline

Then, you can run the pipeline by running:
`nextflow run ITS-pipeline -profile manoa_hpc --reads 'PATH_TO_READS/GLOB_PATTERN'`

To run the pipeline using Docker, you need to first create the docker container:
```
> make python_container
> make R_container
```

Then, you can run the pipeline using the docker profile

## Pipeline summary

The ITS analysis pipeline is summarized below. Values in curly braces ({}) correspond to default values of tunable parameters.

**Inputs**: 
Reads need to be demultiplexed and gzipped

**ITS extraction**: 
Target region ({ITS1}/ITS2) is extracted using ITSxpress. Reads missing either the left or the right flanking regions are discarded. If reads are paired, both reads are merged before ITS region extraction. Since reverse reads are often of very low quality, default mode is single-end.

**Contig filtering (Python, Fastx-toolkit, VSEARCH)**: 
Contigs containing 'N' nucletodes are discarded, as well as contigs smaller than {20} bp. Then, we use fastq_quality_filter and keep reads with at worst {90%} of their bases above quality {25}. We further filter the contigs and remove any chimeric sequence using VSEARCH.

**Denoising (Dada2)**: 
`learnErrors()`, `dada()`: Error models and denoising are performed on each sample independently.

**OTU clustering (VSEARCH)**: 
OTU are clustered at similarity levels {100, 97}% (100% means no clustering).

**Co-occurrence pattern correction (LULU)**: 
A daughter OTU is merged with its parent if:
* they share at least {97}% similarity
* {min}(daughter\_abundance\_sample/parent\_abundance\_sample) < {1}
* the relative co-occurence (proportion of time the daughter is present when the parent is present) must be at least {1}

**Consensus classification (VSEARCH)**: 
Lineages are assigned to each individual sequence using the UNITE reference database. Consensus taxonomy is done for each OTU using a consensus vote. If the consensus is lower than {50%} as a given rank, the taxonomy is not reported.

**Summaries**: 
(samples x pipeline steps) table with the number of remaining sequences in each sample at each step
