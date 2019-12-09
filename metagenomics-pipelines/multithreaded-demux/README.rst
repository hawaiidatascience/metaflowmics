Read demultiplexing
===================

Pre-requisites
--------------

- `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_
- python3 + libraries: Biopython, pandas, numpy, matplotlib, seaborn, argparse

Usage
-----

Clone the repository:

.. code-block:: bash

    git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
	cd nextflow_cmaiki/metagenomics-pipelines

Running the pipeline
^^^^^^^^^^^^^^^^^^^^

To run the pipeline on your data, simply enter the following command:

.. code-block:: bash

    nextflow run multithreaded-demux -profile <config> --inputdir "<path_to_reads>"

For more information about the available profiles, see the corresponding section.

Demultiplexing steps
--------------------

The demultiplexing algorithm is summarized below. Values in between "<" ">" correspond to default values of tunable parameters.

Inputs
^^^^^^
This algorithm is expecting 5 input files for demultiplexing:

- 2 index fastq files (unzipped) matching the glob pattern `*_I{1,2}*.fastq.gz`
- 2 read fastq files (unzipped) matching the glob pattern `*_R{1,2}*.fastq.gz`
- 1 barcode file (extension: .csv), comma separated, with no header and 3 columns: (sample name, forward barcode, reverse complement of reverse barcode)

Guessing the mapping order
^^^^^^^^^^^^^^^^^^^^^^^^^^
To demultiplex, we need to map each index pair in the (I1, I2) files to a barcode pair in the barcode file. Due to protocol variability in sample preparation, the matching order is not always consistent, and for two different experiment, we can have either the "direct" matching order (:math:`I1 \equiv barcode1, I2 \equiv barcode2`) or a "reversed" matching order (:math:`I1 \equiv barcode2, I2 \equiv barcode1`). If the user knows the matching order, it is recommended to enfore it by setting the --matching flag to either "direct" or "reversed". 

The default behavior is a third option, "auto". The guess consists in counting the occurrences of each index pairs from (I1, I2) (in this order), extract the 20 most frequent pairs, and see how many of those eactly match a given barcode pair. The same is done for the pair (I2, I1). If (I1, I2) has the most matches, then the order is set to "direct", and otherwise to "reversed". To make the process faster, any index pair containing a N nucleotides are discarded.

Mapping index to sample names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once the order is determined, both index and sequencing reads are split in chunks of equal sizes of ceil(n_reads/{nsplits}). For each file chunk, we compare each index pair with all barcodes and extract the samples with less errors than {max_mismatches} for each forward and reverse match (therefore in total, there can be up to :math:`2*max_mismatches-1` errors). If several samples match this criteria, we keep the sample with the least mismatches.

Figures
^^^^^^^
To help the user choose parameters for the downstream analysis, two figures are provided:

- The first figure is the FASTQC html plot showing the base quality distribution along the read. This plot is calculated on the pooled samples for the forward and reverse reads respectively.
- The second figure is the log sample size distribution.
