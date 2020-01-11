Microbial 16S pipeline
======================

Pre-requisites
--------------

- `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_
- python(>=3.6) + libraries: Biopython, pandas, matplotlib, seaborn
- R(>=3.5) + libraries: ggplot2, lulu, dada2, seqinr, stringr, ShortRead, doParallel, ape, phyloseq
- `Mothur <https://github.com/mothur/mothur>`_ (tested with v1.43) 

Download the software
^^^^^^^^^^^^^^^^^^^^^

Clone the repository:

.. code-block:: bash

    git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
    cd nextflow_cmaiki/metagenomics-pipelines

Silva Database
^^^^^^^^^^^^^^

In addition, you will need to download the Silva reference database available on the `Mothur website <https://www.mothur.org/wiki/Silva_reference_files>`_:

.. code-block:: bash

	mkdir -p databases && cd databases
	wget https://www.mothur.org/w/images/3/32/Silva.nr_v132.tgz
	tar -xvzf Silva.nr_v132.tgz
	rm -f Silva.nr_v132.tgz && cd ..

For species assignments, the database is already available in the repository at `metagenomics-pipelines/databases/silva_species_assignment_v132.fa.gz`
	
Usage
-----

To run the pipeline on your data, simply enter the following command:

.. code-block:: bash

    nextflow run 16S-pipeline -profile <config> --reads "<path_to_reads/glob_pattern>" --referenceAln databases/silva.nr_v132.align --referenceTax databases/silva.full_v132.tax

The input reads need to be in the `.fastq` format (preferably gzipped) in a single folder. Reads can be single or paired-end. In the former case, the flag `--singleEnd` must be set and in the latter case, the glob pattern needs to group the R1 and R2 reads using the syntax `*R{1,2}*`. 
	
For more information about the available profiles, see the `profiles <https://metagenomics-pipelines.readthedocs.io/en/latest/getting_started.html#configuration-profiles>`_ section.

Outputs
-------

After running the pipeline, the results are stored in the output directory you chose (default is 16S_pipeline_outputs). The results are organized the following way:

- The folder `Results/` groups all of pipeline final outputs:

  - The subfolder `figures/` provides useful insights into the data outputs.
  - The subfolder `main/` contains the pipeline data that is ready for further analysis, and includes:

	- The abundance table: A (sample x OTU) abundance table (tsv formatted) with the `.shared` extension.
	- The taxonomy table: A (OTU x 6 taxonomic tanks) table (tsv formatted) with the `.taxonomy` extension.
	- The OTU sequences: A fasta formatted file with the representative sequences for each OTU.
	- A mothur "database" file that combines all of the previous 3 files.

  - The subfolder `raw/` contains the same information as `main`, before the subsampling, taxa filtering, multipletons filter and lulu steps.
  - The subfolder `postprocessing/` provides alpha and beta diversity metrics. For more details, see `16S postprocessing. <https://metagenomics-pipelines.readthedocs.io/en/latest/pipeline_16S.html#postprocessing>`_

- The folder `Misc/` is mainly for troubleshooting purposes. It groups all of the intermediate results from the pipeline, one sub-folder per step. 

16S pipeline steps
------------------

The 16S analysis pipeline is summarized below. Values in between "<" ">" correspond to default values of tunable parameters.

Inputs
^^^^^^
Reads need to be demultiplexed and gzipped

Read filtering and denoising (Dada2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reads are truncated at positions *<220> / <190>* or at the first occurrence of a base of quality *<2>* or lower. Reads matching the phiX genome are *<discarded>*, as well as reads with an expected number of errrors above *<3>*. Reads shorter than *<20bp>* are filtered out. Finally, samples with less than *<50>* reads are discarded.
Error models and denoising are then performed on each sample independently

Read merging (Dada2)
^^^^^^^^^^^^^^^^^^^^
Paired reads are merged if they overlap by at least *<20bp>* with *<1bp>* mismatch at most.

Contig filtering (Mothur)
^^^^^^^^^^^^^^^^^^^^^^^^^
Contigs are aligned against the silva reference database. Discard any sequence with an alignment shorter than *<50bp>*, as well as sequences starting after where *<95%>* of the sequences start, or end before *<95%>* of the sequences end.

Chimeric contigs are removed using Mothur's implementation of VSEARCH.

OTU clustering (Mothur)
^^^^^^^^^^^^^^^^^^^^^^^
OTU are clustered at similarity levels *<100%, 97%>* (100% means no clustering). 

Taxa filter
^^^^^^^^^^^
Lineages are assigned to each individual sequence using the SILVA reference database. Any sequence matching *<mitochondria, chloroplasts, unknown>* annotations are removed.

Multipletons filter
^^^^^^^^^^^^^^^^^^^
OTU with a total abundance of *<2>* or below are discarded.

.. _subsampling:

Subsampling
^^^^^^^^^^^
We perform sample normalization by subsampling each sample to the same level. Samples with a size below this level are discarded. By default, the subsampling level is defined as the *<10th>* percentile of the sample sizes, and a hard threshold is set if this value goes below *<5000>*. The recommended approach is to determine this value before the analysis and a custom subsampling level can be set. This step can be skipped.

Co-occurrence pattern correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A daughter OTU is merged with its parent if:

* they share at least *<97%>* similarity
* `daughter_abundance < parent_abundance`: in all samples (*<"min">*) or in average ("avg").
* the relative co-occurence (proportion of time the daughter is present when the parent is present) must be at least *<1>*

Rare sequences filter
^^^^^^^^^^^^^^^^^^^^^
OTU with a total abundance of *<2>* or below are discarded.

Consensus classification and representative sequences extraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using the remaining sequences, we choose a representative sequence for each OTU cluster as the most abundant sequence in the cluster. 
For each taxonomic rank, OTU's taxonomy is assigned as the majority vote in the OTU cluster. If the consensus vote is lower than 51%, no taxonomy is assigned at the given rank.

Summaries
^^^^^^^^^
- (samples x pipeline steps) table with the number of remaining sequences in each sample at each step
- Figures

  #. (top OTUs x samples) bi-clustered heatmap with phylum, class and order information.
  #. scatter plot of OTUs abundance vs prevalence, one facet per phylum.
  #. scatter plot of OTUs abundance vs prevalence for proteobacteria, one facet per class.
  #. barplot of relative taxonomy composition at Phylum level for each sample. In a metadata table is provided, this plots represents the composition for each level of the provided factor.

.. _16Spostproc:
	 
Postprocessing
^^^^^^^^^^^^^^
For each clustering threshold, we compute alpha and beta diversity metrics (see `mothur calculators <https://www.mothur.org/wiki/Calculators>`_ for a full description of these acronyms)

- Alpha diversity: `nseqs`, `sobs`, `chao`, `shannon`, `shannoneven`
- Beta diversity: `braycurtis`, `thetayc`, `sharedsobs`, `sharedchao`

In addition, we compute the phylogenetic tree using `FastTree <http://www.microbesonline.org/fasttree/>`_ and compute the UniFrac distances using the R's `phyloseq <https://bioconductor.org/packages/release/bioc/html/phyloseq.html>`_ package implementing the `Fast UniFrac <https://www.ncbi.nlm.nih.gov/pubmed/19710709>`_ algorithm.
