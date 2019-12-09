# 16S-pipeline

## Getting started

Follow these instructions to get the pipeline started on your machine

### Pre-requisites

The preferred way to use this pipeline is through the singularity configuration. It requires little setup since all dependencies are dealt with using a custom docker container.

Requirements:
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

If you wish to run the pipeline without it, you will need to satisfy in addition the following dependencies:

- python3 + libraries: Biopython, pandas, matplotlib, seaborn
- R(>=3.5) + libraries: ggplot2, lulu, dada2, seqinr, stringr, ShortRead, doParallel, ape, phyloseq
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

You can also run the pipeline using the available docker instances. For that, you'll just need to select the correct profile "docker"
```
nextflow run 16S-pipeline -profile docker --reads 'PATH_TO_READS/GLOB_PATTERN'
```

Then, you can run the pipeline using the docker profile

## Pipeline summary

The 16S analysis pipeline is summarized below. Values in curly braces ({}) correspond to default values of tunable parameters.

**Inputs**: 
Reads need to be demultiplexed and gzipped

**Read filtering (Dada2)**: 
`filterAndTrim()`: Reads are truncated at positions {220} / {190} (fwd/rev) or at the first occurrence of a base of quality {2} or lower. Reads matching the phiX genome are {discarded}, as well as reads with an expected number of errrors above {maxEE}. Reads shorter than {20} bp are filtered out. Finally, samples with less than 50 reads are discarded.

**Denoising (Dada2)**: 
`learnErrors()`, `dada()`: Error models and denoising are performed on each sample independently.

**Read merging (Dada2)**: 
`mergePairs()`: Paired reads are merged if they overlap by at least {20} bp with {1} mismatch at most

**Contig filtering (Mothur)**: 
Contigs are aligned against the silva reference database. Discard any sequence with an alignment shorter than {50} bp, as well as sequences starting after where {95}% of the sequences start, or end before {95}% of the sequences end.

**Chimera filtering (Mothur / VSEARCH)**: 
Chimeric contigs are removed using Mothur's implementation of VSEARCH

**OTU clustering (Mothur)**: 
OTU are clustered at similarity levels {100, 97}% (100% means no clustering). 

**Taxa filter**: 
Lineages are assigned to each individual sequence using the SILVA reference database. Any sequence matching {mitochondria, chloroplasts, unknown} annotations are removed.

**Multipletons filter**: 
OTU with a total abundance of {2} or below are discarded.

**Subsampling**: 
We perform sample normalization by subsampling each sample to the same level. Samples with a size below this level are discarded. By default, the subsampling level is defined as the {10th} percentile of the sample sizes, and a hard threshold is set if this value goes below {5000}. The recommended approach is to determine this value before the analysis and a custom subsampling level can be set. This step can be skipped.

**Co-occurrence pattern correction**: 
A daughter OTU is merged with its parent if:
* they share at least {97}% similarity
* {min}(daughter\_abundance\_sample/parent\_abundance\_sample) < {1}
* the relative co-occurence (proportion of time the daughter is present when the parent is present) must be at least {1}

**Rare sequences filter**: 
OTU with a total abundance of {2} or below are discarded.

**Consensus classification and representative sequences extraction**
Using the remaining sequences, we choose a representative sequence for each OTU cluster as the most abundant sequence in the cluster. 
For each taxonomic rank, OTU's taxonomy is assigned as the majority vote in the OTU cluster. If the consensus vote is lower than 51%, no taxonomy is assigned at the given rank.

**Summaries**: 
- (samples x pipeline steps) table with the number of remaining sequences in each sample at each step
- Figures
  - (top OTUs x samples) bi-clustered heatmap with phylum, class and order information.
  - scatter plot of OTUs abundance vs prevalence, one facet per phylum.
  - scatter plot of OTUs abundance vs prevalence for proteobacteria, one facet per class.
  - barplot of relative taxonomy composition at Phylum level for each sample. In a metadata table is provided, this plots represents the composition for each level of the provided factor.

**Postprocessing**: 
For each clustering thresho, we compute alpha and beta diversity metrics (see [mothur calculators](https://www.mothur.org/wiki/Calculators) for a full description of these acronyms)
- Alpha diversity: `nseqs`, `sobs`, `chao`, `shannon`, `shannoneven`
- Beta diversity: `braycurtis`, `thetayc`, `sharedsobs`, `sharedchao`
In addition, we compute the phylogenetic tree using [FastTree](http://www.microbesonline.org/fasttree/) and compute the UniFrac distances using the R's [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html) package implementing the [Fast UniFrac](https://www.ncbi.nlm.nih.gov/pubmed/19710709) algorithm.
