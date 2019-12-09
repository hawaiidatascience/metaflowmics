Fungal ITS pipeline
===================

Pre-requisites
--------------

- `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_
- python3 + packages: Biopython, pandas, matplotlib, seaborn, ITSxpress
- R(>=3.5) + libraries: ggplot2, lulu, dada2, seqinr, stringr, ShortRead, doParallel, ape, phyloseq
- install `VSEARCH <https://github.com/torognes/vsearch/releases>`_, `HMMER <http://eddylab.org/software/hmmer>`_ and `BBTools <https://sourceforge.net/projects/bbmap>`_

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
    nextflow run ITS-pipeline -profile CONFIG --reads "PATH_TO_READS/GLOB_PATTERN"

For more information about the available profiles, see the corresponding section.

ITS pipeline steps
------------------

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
