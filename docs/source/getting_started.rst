Getting Started
===============

Downloading the pipelines
-------------------------

To download the pipelines, simply clone the repository:

.. code-block:: bash

    git clone https://github.com/hawaiidatascience/nextflow_cmaiki.git
    cd nextflow_cmaiki/metagenomics-pipelines


Docker / Singularity configuration (recommended)
------------------------------------------------

Each of the available pipelines are were built using multiple computing languages and softwares. The dependencies are specified in the appropriate sections. Since setting up these environment can be time time consuming, and sometime complicated depending on the platform, we provide a unified docker instance for all pipelines. This instance can be then used by `Docker` or `Singularity` (if the user does not have sudo rights). This is th recommened way to use the pipeline.

The requirements to run the pipeline with the docker instance are minimal:

- `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_ (>=19.10)
- Singularity can be compiled from source in the `GitHub releases <https://github.com/sylabs/singularity/releases>`_ or installed via `apt-get` in the NeuroDebian repository.
- For the 16S and ITS pipelines, downloading the microbial / fungal database is also required. Refer to the corresponding sections for more information.


Usage
-----

The pipelines have a similar functionning. To get the available options, you can display the help using the following command (where <pipeline-name> refers to one of the 3 folders in `nextflow_cmaiki/metagenomics-pipelines`):

.. code-block:: bash

    nextflow run <pipeline-name> --help

When running any of the pipelines, you will need to provide a few mandatory arguments:

#. The path to the input reads (see corresponding section for more details)
#. The profile defined in conf/*.config. Several profiles are already available, such as:

   - The `singularity` and `docker` profile (recommended)
   - The `hpc` profile is set up for the slurm scheduler and uses the singularity image. This profile requires `singularity` to be installed on the HPCC as a module. See the advanced configuration section for custom setup (module definition and scheduler information).
   - The `gcp` profile to run the pipelines on the Google Cloud Platform. See the corresponding section for how to set it up.
   - The `local` profile does not use the docker instance and therefore requires the machine to have all the required packages and softwares.

To select a given configuration, you just need to append `-profile <profile_name>` to the command.
For custom configuration, please see the advanced configuration section.
