Documentation: Nextflow pipelines for metagenomics
==================================================
.. image:: https://travis-ci.org/hawaiidatascience/nextflow_cmaiki.svg?branch=master
    :target: https://travis-ci.org/hawaiidatascience/nextflow_cmaiki

This repository is a collection of tools for 16S and ITS data analysis and was developed to support scientifics at the `Pacific Biosciences Research Center <http://www.pbrc.hawaii.edu/>`_ to analyze microbial data.
This work was funded by the `Hawaii Data Science Institute <http://datascience.hawaii.edu/>`_

Each analysis tool has been implemented using `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_ to facilitate their use on any environment, and more specifically on a High Performance Computing Cluster or a Cloud service (GCP, AWS, ...). Some configurations, including the one for the HPCC at the University of Hawaii at Manoa is already available in the `conf` folder. These configurations can be expanded for any other platform. 
See `nextflow documentation <https://www.nextflow.io/docs/latest/config.html>`_ for more information. 

You will find three available app:

#. A demultiplexing app that can demultiplex reads with single or paired barcodes.
#. A 16S pipeline for bacterial reads
#. An ITS pipeline for fungal reads

If you are interested in using one of these app, please see the README in the respective folder.

Finally, a docker instance is already available. To use them, you can use the "docker" (or "singularity" if you do not have sudo rights) configuration file.

