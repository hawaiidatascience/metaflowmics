Fungal ITS pipeline
===================

Pre-requisites
--------------

- `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_
- python(>=3.6) + packages: Biopython, pandas, matplotlib, seaborn, ITSxpress
- R(>=3.5) + libraries: ggplot2, lulu, dada2, seqinr, stringr, ShortRead, doParallel, ape, phyloseq
- install `VSEARCH <https://github.com/torognes/vsearch/releases>`_, `HMMER <http://eddylab.org/software/hmmer>`_ and `BBTools <https://sourceforge.net/projects/bbmap>`_

Download the software
^^^^^^^^^^^^^^^^^^^^^

Clone the repository:

.. code-block:: bash

    git clone https://github.com/hawaiidatascience/metaflowmics.git
    cd metaflowmics/metaflowmics

Usage
-----

To run the pipeline on your data, simply enter the following command:

.. code-block:: bash

    nextflow run ITS-pipeline -profile <config> --reads "<path_to_reads/glob_pattern>"

The input reads need to be in the `.fastq` format (preferably gzipped) in a single folder. Reads can be single or paired-end. In the latter case, the glob pattern needs to group the R1 and R2 reads using the syntax "\*R{1,2}\*", and the flag `--paired_end` must be set.
	
For more information about the available profiles, see the `profiles <https://metagenomics-pipelines.readthedocs.io/en/latest/getting_started.html#configuration-profiles>`_ section.

ITS pipeline steps
------------------

The ITS analysis pipeline is summarized below. Values in between "<" ">" correspond to default values of tunable parameters.

Inputs
^^^^^^
Reads need to be demultiplexed and gzipped

ITS extraction
^^^^^^^^^^^^^^
The target region *<ITS1>* is extracted using ITSxpress. Reads missing either the left or the right flanking regions are discarded. If reads are paired, both reads are merged before ITS region extraction. Since reverse reads are often of very low quality, default mode is single-end.

Contig filtering (DADA2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Contigs containing 'N' nucletodes are discarded, as well as contigs smaller than *<20bp>*. Then, we use `filterAndTrim()` and keep reads with at most *<3>* expected errors. We further filter the contigs and remove any chimeric sequence using VSEARCH.

Denoising (DADA2)
^^^^^^^^^^^^^^^^^
`learnErrors(), dada()`: Error models and denoising are performed on each sample independently.

OTU clustering (VSEARCH)
^^^^^^^^^^^^^^^^^^^^^^^^
OTU are clustered at similarity levels *<100%, 97%>* (100% means no clustering).

Co-occurrence pattern correction (LULU)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A daughter OTU is merged with its parent if:

* they share at least *<97%>* similarity
* daughter_abundance < parent_abundance in all samples (*<min>*) or in average (*avg*).
* the relative co-occurence (proportion of time the daughter is present when the parent is present) must be at least *<1>*

Classification (DADA2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Lineages are assigned to each individual sequence using the UNITE reference database. The assignment is performed with the naive Bayesian classifier method and sequences with a confidence at a given taxonomic rank lower than <50%> are not assigned.

Summaries
^^^^^^^^^
(samples x pipeline steps) table with the number of remaining sequences in each sample at each step
