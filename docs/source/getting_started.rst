.. _getting_started:

Getting Started
===============

Downloading the pipelines
-------------------------

To download the pipelines, simply clone the repository:

.. code-block:: bash

    git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
    cd nextflow_cmaiki/metagenomics-pipelines

Usage
-----

Each pipeline work the same way. To get the available options, you can display the help using the following command (where <pipeline-name> refers to one of the 3 folders in `nextflow_cmaiki/metagenomics-pipelines`):

.. code-block:: bash

    nextflow run <pipeline-name> --help

When running any of the pipelines, you will need to provide a few mandatory arguments:

#. The path to the input reads (see corresponding section for more details)
#. The profile defined in conf/*.config.

Dependencies
------------

Docker / Singularity configuration (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each of the available pipelines are were built using multiple computing languages and softwares. The dependencies are specified in the appropriate sections. Since setting up these environment can be time time consuming, and sometime complicated depending on the platform, we provide a unified docker instance for all pipelines. This instance can be then used by `Docker` or `Singularity` (if the user does not have sudo rights). This is th recommened way to use the pipeline.

The requirements to run the pipeline with the docker instance are minimal:

- `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_ (>=19.10)
- Singularity can be compiled from source in the `GitHub releases <https://github.com/sylabs/singularity/releases>`_ or installed via `apt-get` in the NeuroDebian repository.
- For the 16S and ITS pipelines, downloading the microbial / fungal database is also required. Refer to the corresponding sections for more information.

Custom configuration
^^^^^^^^^^^^^^^^^^^^

Alternatively, you can install all required dependencies on your machine. The specific requirements for each pipeline is specified on the corresponding documentation. 

Configuration profiles
----------------------

To facilitate the configuration of the pipeline, run profiles are available in the `conf` folder for multiple use cases:

- The `singularity` and `docker` profile to run the pipeline locally
- The `local` profile does not use the docker instance and therefore requires the machine to have all the required packages and softwares.
- The `hpc` profile is set up for the slurm scheduler and uses the singularity image. This profile requires `singularity` to be installed on the HPCC as a module. See the :ref:`hpc_conf` section for more information.
- The `gcp` profile to run the pipelines on the Google Cloud Platform. See the :ref:`gcp_conf` section for more information.

To select a given configuration, you just need to append `-profile <profile_name>` to the command.
For custom configuration, see the :ref:`advanced_config` section.
