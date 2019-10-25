# Nextflow Demultiplexing 

## Getting started

Follow these instructions to get the app started on your machine

### Pre-requisites

To run the app, you will need to satisfy the following dependencies:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- python3 + libraries: Biopython, pandas, numpy, matplotlib, seaborn, argparse

### Usage

- Clone the repository:
```
git clone https://github.com/Puumanamana/nextflow-demux.git
```

#### Make your own configuration file

To run the pipeline locally, you need to set up a configuration file. An example is available in `conf/local.config`.
You can modify this configuration to fit the specs of your machine.

#### Changing the default pipeline parameters

Pipeline parameters can be set either:
- in the file `nextflow.config`
- by using the flags --[PARAMETER-NAME] [PARAMETER-VALUE] (see `nextflow run nextflow-demux --help`)

### Demultiplexing

The demultiplexing algorithm is summarized below. Values in curly braces ({}) correspond to default values of tunable parameters.

**Inputs**: 
This algorithm is expecting 5 input files for demultiplexing:
- 2 index fastq files (unzipped) matching the glob pattern "\*_I{1,2}\*.fastq.gz"
- 2 read fastq files (unzipped) matching the glob pattern "\*_R{1,2}\*.fastq.gz"
- 1 barcode file named "barcodes.csv", comma separated, with no header and 3 columns: 
   (sample name, forward barcode, reverse complement of reverse barcode)

**Guessing the mapping order**: 
To demultiplex, we need to map each index pair in the (I1, I2) files to a barcode pair in the barcode file. Due to protocol variability in sample preparation, the matching order is not always consistent, and for two different experiment, we can have either the "direct" matching order (I1 <-> barcode1, I2 <-> barcode2) or a "reversed" matching order (I1 <-> barcode2, I2 <-> barcode1). If the user knows the matching order, it is recommended to enfore it by setting the --matching flag to either "direct" or "reversed". 

The default behavior is a third option, "auto". The guess consists in counting the occurrences of each index pairs from (I1, I2) (in this order), extract the 20 most frequent pairs, and see how many of those eactly match a given barcode pair. The same is done for the pair (I2, I1). If (I1, I2) has the most matches, then the order is set to "direct", and otherwise to "reversed". To make the process faster, any index pair containing a N nucleotides are discarded.

**Mapping index to sample names**: 
Once the order is determined, both index and sequencing reads are split in chunks of equal sizes of ceil(n_reads/{nsplits}). For each file chunk, we compare each index pair with all barcodes and extract the samples with less errors than {max_mismatches} for each forward and reverse match (therefore in total, there can be up to 2*{max_mismatches}-1 errors). If several samples match this criteria, we keep the sample with the least mismatches. If there are still multiple matches, we discard the read pair. Other strategy are available to handle multi-mappers:
- The "discard" strategy discards any read pair that maps multiple samples. Depending on the {max_mismatches} value, this can result in most of the reads being discarded.
- The "min_all" strategy keep all matches achieving the minimum score. Therefore, if a read pair achieves the minimum score for multiple samples, it is copied in each.
- The "all" strategy keeps all matches. Therefore it copies the read pair in all the matching samples.

**Figures**: 
To help the user choose parameters for the downstream analysis, two figures are provided:
- The first figure is the FASTQC html plot showing the base quality distribution along the read. This plot is calculated on the pooled samples for the forward and reverse reads respectively.
- The second figure is the log sample size distribution.
