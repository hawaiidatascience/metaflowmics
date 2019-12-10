.. _advanced_config:

Computing resources
===================

We saw that we can change the default parameters of any application by modifying the `nextflow.config` file. However, the computing resources are defined separately in the `conf` folder and they should be tuned to your work environment. You will find there multiple pre-filled configuration to help you. If you are working on a remote server or your local computer, the recommended configuration is the `conf/docker` configuration or, if you do not have sudo access, the `conf/singularity` configuration. Otherwise, configuration for a HPCC or Google Cloud Platform are provided.

Main parameters
---------------

- queueSize
- cpus
- memory
- time
- errorStrategy
- maxRetries
  
.. include:: hpcc.rst

.. include:: gcp.rst
			 
