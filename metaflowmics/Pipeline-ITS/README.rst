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

    git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
    cd nextflow_cmaiki/metagenomics-pipelines

UNITE Database
^^^^^^^^^^^^^^

In addition, you will need to download the UNITE reference database (all eukaryotes) available on the `UNITE website <https://unite.ut.ee/repository.php>`_. Converting any non ASCII character is also preferable using the `iconv` command:

.. code-block:: bash

    mkdir -p databases && cd databases
    wget https://files.plutof.ut.ee/doi/A1/C9/A1C964DFB03C2A1B37FA16784BA739C88F0941AC68560CEA54DD707F1CF00AC4.zip -O uniteDB.zip
    unzip uniteDB.zip && rm uniteDB.zip
    iconv -f utf-8 -t ascii//translit UNITE_public_all_02.02.2019.fasta > uniteDB.fasta
    rm -f UNITE_public_all_02.02.2019.fasta && cd ..

Usage
-----

To run the pipeline on your data, simply enter the following command:

.. code-block:: bash

    nextflow run ITS-pipeline -profile <config> --reads "<path_to_reads/glob_pattern>" --uniteDB databases/uniteDB.fasta

The input reads need to be in the `.fastq` format (preferably gzipped) in a single folder. Reads can be single or paired-end. In the latter case, the glob pattern needs to group the R1 and R2 reads using the syntax "\*R{1,2}\*", and the flag `--pairedEnd` must be set.
	
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

Contig filtering (Python, Fastx-toolkit, VSEARCH)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Contigs containing 'N' nucletodes are discarded, as well as contigs smaller than *<20bp>*. Then, we use fastq_quality_filter and keep reads with at worst *<90%>* of their bases above quality *<25>*. We further filter the contigs and remove any chimeric sequence using VSEARCH.

Denoising (Dada2)
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

Consensus classification (VSEARCH)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Lineages are assigned to each individual sequence using the UNITE reference database. Consensus taxonomy is done for each OTU using a consensus vote. If the consensus is lower than *<50%>* as a given rank, the taxonomy is not reported.

Summaries
^^^^^^^^^
(samples x pipeline steps) table with the number of remaining sequences in each sample at each step
