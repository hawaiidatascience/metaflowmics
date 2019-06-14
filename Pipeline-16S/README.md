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
```
git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
cd Pipeline-16S
```

#### Make your own configuration file
To run the pipeline locally, you need to set up a configuration file. An example is available in `conf/poire.config`.
Once the configuration file is set, you can add the entry in the `nextflow.config file` the following way:
Under "profiles", add the lines
``` 
[MY_CONFIG_NAME]{
    includeConfig 'conf/[PATH_TO_MY_CONFIG]'
}
```

#### Changing the default pipeline parameters

Pipeline parameters can be set either:
- in the file `nextflow.config`
- by using the flags --[PARAMETER-NAME] [PARAMETER-VALUE] (see `nextflow run 16S-pipeline --help`)

### Running the pipeline

Then, you can run the pipeline by running:
`nextflow run 16S-pipeline -profile manoa_hpc --reads 'PATH_TO_READS/GLOB_PATTERN'`

## Pipeline summary

The 16S analysis pipeline is summarized below and includes, for each step, the tunable parameters and their default values.

**dada2**
- filterAndTrim > Parameters: 
  - minimum read length=20
  - maximum expected errors=3
  - read length truncation=(fwd: 220, rev: 190)
- learnErrors
- mergePairs > Parameters: 
  - minimum overlap=20
  - maximum mismatches=1

**Mothur**
--> alignment against reference for further filtering
- align.seqs; filter.seqs; screen.seqs > Parameters: 
  - criteria=95 (remove any sequence that starts after the position that [CRITERIA]% of the sequences do, or ends before the position that [CRITERIA]% of the sequences do)
- chimera.vsearch; remove.seqs
- sub.sample (if skipSubSampling=false) > Parameters: 
  - subsamplingQuantile=0.10 (subsample at the 10th percentile of the sample sizes)
  - minSubsampling=5000 (subsample at 5k if the 10th percentile is below 5k)
- cluster (at 95,97,99 and 100% identity)
- classify.seqs ; classify.otu 

**Lulu** (see https://rdrr.io/github/tobiasgf/lulu/man/lulu.html) for more details)
- preLulu (create matchlists for LULU)
- LULU > Parameters: 
  - min\_ratio\_type="min" (function to aggregate the daughter/parent abundance ratios across samples where the parent OTU is present)
  - min\_ratio=1 (a daughter is merged with the parent if the daughter is at most [min\_ratio] times less abundant than the parent)
  - min\_match=97
  - min\_rel\_cooccurence=1 (the daughter must be present each time the parent is present)

In summary, a daughter OTU is merged with its parent if:
* they share at least [min\_match]% similarity
* min(daughter\_abundance\_sample/parent\_abundance\_sample) < [min\_ratio]
* the relative co-occurence (proportion of time the daughter is present when the parent is present) must be at least [min\_rel\_cooccurence] 

**Python**
- OutputFilter (remove LULU flagged sequences from FASTA and taxonomy file; remove low abundant OTUs (<=2))

**Postprocessing**
- ConvertToMothur: Convert output file to .shared file
- Results: Mothur postprocessing with get.relabund, clearcut and unifrac.weighted
- SummaryFile: Generates summary of reads per sample per step
