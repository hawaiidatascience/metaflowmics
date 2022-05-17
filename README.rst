MetaFlow|mics: Nextflow pipelines for metagenomics
==============================================
.. image:: https://readthedocs.org/projects/metagenomics-pipelines/badge/?version=latest
   :target: https://metagenomics-pipelines.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/github/workflow/status/hawaiidatascience/metaflowmics/test?label=tests&color=orange&logo=github&logoColor=orange
   :target: https://github.com/hawaiidatascience/metaflowmics/actions?query=workflow

Cite
----

Cédric Arisdakessian, Sean B. Cleveland, and Mahdi Belcaid. 2020. MetaFlow|mics: Scalable and Reproducible Nextflow Pipelines forthe Analysis of Microbiome Marker Data. InPractice and Experience in Advanced Research Computing (PEARC ’20), July 26–30, 2020,Portland, OR, USA.ACM, New York, NY, USA, 9 pages.

Description
-----------

This repository is a collection of tools for 16S and ITS data analysis and was developed to support scientifics at the `Pacific Biosciences Research Center <http://www.pbrc.hawaii.edu/>`_ to analyze microbial data.
This work was funded by the `Hawaii Data Science Institute <http://datascience.hawaii.edu/>`_.

Each analysis tool has been implemented using `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_ to facilitate their use on any environment, and more specifically on a High Performance Computing Cluster or a Cloud service (GCP, AWS, ...). Some configurations are already available in the `conf` folder (local, hpcc, google cloud). These configurations can be expanded for any other platform. 

You will find three available app:

#. A demultiplexing app for single or paired barcodes.
#. A 16S pipeline for bacterial reads
#. An ITS pipeline for fungal reads

See the `documentation <https://metagenomics-pipelines.readthedocs.io>`_ for more details.

Contribute
----------
- Issue tracker: `GitHub <https://github.com/hawaiidatascience/metaflowmics/issues>`_
- Source Code: `GitHub <https://github.com/hawaiidatascience/metaflowmics/tree/master/metaflowmics>`_

License
-------
Apache 2.0
