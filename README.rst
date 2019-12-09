Nextflow pipelines for metagenomics
===================================
.. image:: https://travis-ci.org/hawaiidatascience/nextflow_cmaiki.svg?branch=master
   :target: https://travis-ci.org/hawaiidatascience/nextflow_cmaiki
.. image:: https://readthedocs.org/projects/metagenomics-pipelines/badge/?version=latest
   :target: https://metagenomics-pipelines.readthedocs.io/en/latest/?badge=latest
			 
Cite (Work in progress)
-----------------------

Scalable and reproducible pipelines for the analysis of microbiome marker data

Description
-----------

This repository is a collection of tools for 16S and ITS data analysis and was developed to support scientifics at the `Pacific Biosciences Research Center <http://www.pbrc.hawaii.edu/>`_ to analyze microbial data.
This work was funded by the `Hawaii Data Science Institute <http://datascience.hawaii.edu/>`_.

Each analysis tool has been implemented using `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_ to facilitate their use on any environment, and more specifically on a High Performance Computing Cluster or a Cloud service (GCP, AWS, ...). Some configurations are already available in the `conf` folder (local, hpcc, google cloud). These configurations can be expanded for any other platform. 

You will find three available app:

#. A demultiplexing app for single or paired barcodes.
#. A 16S pipeline for bacterial reads
#. An ITS pipeline for fungal reads

See the `documentation <https://metagenomics-pipelines.readthedocs.io>` for more details.

Contribute
----------
- Issue tracker: `GitHub <https://github.com/hawaiidatascience/nextflow_cmaiki/issues>`
- Source Code: `GitHub <https://github.com/hawaiidatascience/nextflow_cmaiki>`

License
-------
Apache 2.0
