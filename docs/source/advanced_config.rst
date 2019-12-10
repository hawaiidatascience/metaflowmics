.. _advanced_config:

Computing resources
===================

We saw that we can change the default parameters of any application by modifying the `nextflow.config` file. However, the computing resources are defined separately in the `conf` folder and they should be tuned to your work environment. You will find there multiple pre-filled configuration to help you. If you are working on a remote server or your local computer, the recommended configuration is the `conf/docker` configuration or, if you do not have sudo access, the `conf/singularity` configuration. Otherwise, configuration for a HPCC or Google Cloud Platform are provided.

Main parameters
---------------

- queueSize: maximum number of jobs to launch in parallel
- `errorStrategy <https://www.nextflow.io/docs/latest/process.html#errorstrategy>`_: Error handling. Available choices are `terminate`, `finish`, `ignore` and `retry`. Note that some errors (such as segfaults) can sometime be hard to recover from, and the pipeline might exit abruptly even with this parameter set to ignore.
- maxRetries: When `errorStrategy` is set to retry, it corresponds to the maximum number of times a process instance can be re-submitted in case of failure.

In addition, the pipelines define three categories of processes: low / medium / high computation. These resources can be set for either type of resource by changing their values in the corresponding fields. For example, if you wish to change the number of cpus for high computation processeses, you just need to change:

.. code-block:: groovy

   withLabel: high_computation {
      cpus = 3 // Change the value here
   }

- cpus: Number of cpus
- memory: Maximum amount of memory
- time: Maximum running time

.. include:: hpcc.rst

.. include:: gcp.rst
			 
