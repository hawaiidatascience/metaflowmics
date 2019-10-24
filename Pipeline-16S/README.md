# 16S-pipeline

## Getting started

Follow these instructions to get the pipeline started on your machine

### Pre-requisites

To run the pipeline, you will need to satisfy the following dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- python3 + libraries: Biopython, pandas, matplotlib, seaborn
- R + libraries: ggplot2, lulu, dada2, seqinr, stringr, ShortRead
- [Mothur](https://github.com/mothur/mothur) (tested with v1.43) 

### Usage

- Clone the repository:
```
git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
cd Pipeline-16S
```

#### Make your own configuration file
To run the pipeline locally, you need to set up a configuration file. An example is available in `conf/local.config`.
You can modify this configuration to fit the specs of your machine.

#### Changing the default pipeline parameters

Pipeline parameters can be set either:
- in the file `nextflow.config`
- by using the flags --[PARAMETER-NAME] [PARAMETER-VALUE] (see `nextflow run 16S-pipeline --help`)

### Running the pipeline

Then, you can run the pipeline by running:
`nextflow run 16S-pipeline -profile manoa_hpc --reads 'PATH_TO_READS/GLOB_PATTERN'`

To run the pipeline using Docker, you need to first create the docker container:
```
> make python_container
> make R_container
> mothur_container
```

Then, you can run the pipeline using the docker profile

## Pipeline summary

The 16S analysis pipeline is summarized below. Values in curly braces ({}) correspond to default values of tunable parameters.

**Inputs**
- Reads need to be demultiplexed and gzipped

**Read filtering (Dada2)**
`filterAndTrim()`: Reads are truncated at positions {220} / {190} (fwd/rev) or at the first occurrence of a base of quality {2} or lower. Reads matching the phiX genome are {discarded}, as well as reads with an expected number of errrors above {maxEE}. Reads shorter than {20} bp are filtered out. Finally, samples with less than 50 reads are discarded.

**Denoising (Dada2)**
`learnErrors()`, `dada()`: Error models and denoising are performed on each sample independently.

**Read merging (Dada2)**
`mergePairs()`: Paired reads are merged if they overlap by at least {20} bp with {1} mismatch at most

**Contig filtering (Mothur)**
Contigs are aligned against the silva reference database. Discard any sequence with an alignment shorter than {50} bp, as well as sequences starting after where {95}% of the sequences start, or end before {95}% of the sequences end.

**Chimera filtering (Mothur / VSEARCH)**
Chimeric contigs are removed using Mothur's implementation of VSEARCH

**OTU clustering (Mothur)**
OTU are clustered at similarity levels {100, 97}% (100% means no clustering). 

**Consensus classification and taxa filter**
Lineages are assigned to each individual sequence using the SILVA reference database. Consensus taxonomy is done for each OTU and taxa matching {mitochondria, chloroplasts, unknown} are removed.

**Multipletons filter**
OTU with a total abundance of {2} or below are discarded.

**Subsampling**
We perform sample normalization by subsampling each sample to the same level. Samples with a size below this level are discarded. By default, the subsampling level is defined as the {10th} percentile of the sample sizes, and a hard threshold is set if this value goes below {5000}. The recommended approach is to determine this value before the analysis and a custom subsampling level can be set. This step can be skipped.

**Co-occurrence pattern correction**
A daughter OTU is merged with its parent if:
* they share at least {97}% similarity
* {min}(daughter\_abundance\_sample/parent\_abundance\_sample) < {1}
* the relative co-occurence (proportion of time the daughter is present when the parent is present) must be at least {1}

**Rare sequences filter**
OTU with a total abundance of {2} or below are discarded.

**Summaries**
- (samples x pipeline steps) table with the number of remaining sequences in each sample at each step
- Figures

**Postprocessing**
For each clustering thresho, we compute alpha and beta diversity metrics (see [mothur calculators](https://www.mothur.org/wiki/Calculators) for a full description of these acronyms)
- Alpha diversity: `nseqs`, `sobs`, `chao`, `shannon`, `shannoneven`
- Beta diversity: `braycurtis`, `thetayc`, `sharedsobs`, `sharedchao`
