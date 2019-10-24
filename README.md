# App for 16S and ITS data processing

This repository is a collection of tools for 16S and ITS data analysis and was developed to support scientifics at the [Pacific Biosciences Research Center](http://www.pbrc.hawaii.edu/) to analyze microbial data.
This work was funded by the [Hawaii Data Science Institute](http://datascience.hawaii.edu/)

Each analysis tool has been implemented using [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) to facilitate their use on any environment, and more specifically on a High Performance Computing Cluster. Some configurations, including the one for the HPCC at the University of Hawaii at Manoa is already available in the `conf` folder. These configurations can be expanded for any other platform. See [nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information. 

You will find three available apps:
- A demultiplexing app that can demultiplex paired-end reads with paired-barcodes.
- Two pipelines for marker gene sequences processing:
  - A 16S pipeline for bacterial reads
  - An ITS pipeline for fungal reads

If you are interested in using one of these app, please see the README in the respective folder.

Finally, some docker instances are already available with the required environments. To use them, you can use the "docker" configuration file.

