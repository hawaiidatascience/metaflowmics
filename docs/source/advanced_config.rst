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

.. _hpc_conf:

Configure for a HPCC
--------------------

HPCC require a few more specific parameters to be set:

- `executor`: The HPCC scheduler. Default is SLURM, but many others are available. For a comprehensive list see `nextflow documentation <https://www.nextflow.io/docs/latest/executor.html>`_
- `queue`: List of queues where you wish to submit the processes. Multiple queues are separated with commas.
- `module`: Singularity module location.

.. note::

   If you are working on a squashed root system, you will need to pull the repository before running the pipeline. This is likely the case if you get the error `Failed to pull singularity image` when trying to run the pipeline. You will also need to set up a docker account.

   .. code-block:: bash

	  module load path/to/singularity/module
	  mkdir -p ~/.singularity_images.cache
	  mkdir -p /tmp/sg
	  singularity pull --docker-login /tmp/sg/nakor-pipeline-env.img docker://nakor/pipeline-env
	  singularity pull --docker-login /tmp/sg/alexcoppe-fastx.img docker://alexcoppe/fastx
	  singularity pull --docker-login /tmp/sg/pegi3s-fasttree.img docker://pegi3s/fasttree
	  mv /tmp/sg/*.img ~/.singularity_images.cache

.. _gcp_conf:

Configure for Google Cloud Plateform (GCP)
------------------------------------------

In order to use Google Cloud Platform, you need to set up an account on the `Google Cloud website <https://console.cloud.google.com/>`_. Then, you will need to:

#. Create a `project <https://cloud.google.com/resource-manager/docs/creating-managing-projects>`_
#. Create a `storage bucket <https://cloud.google.com/storage/docs/creating-buckets>`_

The last thing you need is to provide nextflow a way to access your account via a token. The instruction for creating the token are extracted from `nextflow documentation <https://www.nextflow.io/docs/latest/google.html>`_:

- Open the `Google Cloud Console <https://console.cloud.google.com/>`_
- Go to APIs & Services → Credentials
- Click on the Create credentials (blue) drop-down and choose Service account key, in the following page
- Select an existing Service account or create a new one if needed
- Select JSON as Key type
- Click the Create button and download the JSON file giving a name of your choice e.g. `creds.json`

Finally define the following variable replacing the path in the example with the one of your credentials file just downloaded:

.. code-block:: bash

   export NXF_MODE=google
   export GOOGLE_APPLICATION_CREDENTIALS=/path/your/file/creds.json

You are almost set to use the GCP platform. The last remaining thing you need to do is set the configuration file in `conf/gcp.conf` and set your project information, the location and type of machines you want to use:

- `project`: name of your project on GCP
- `zone`: Location of the cloud instances. See `GCP documentation <https://cloud.google.com/compute/docs/regions-zones/?hl=en>`_ for a list of available locations.
- `machineType`: specs of the instance you wish to use. A list of available instances is available on the `GCP documentation <https://cloud.google.com/compute/all-pricing>`_.

For more details about the GCP configuration, see the `nextflow documentation <https://www.nextflow.io/docs/latest/google.html>`_

Configure for Amazon Web Service
--------------------------------

Not implemented yet. If you are interested to set it, have a look at the `nextflow documentation <https://www.nextflow.io/docs/latest/awscloud.html>`_
			 
