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
   
   export GOOGLE_APPLICATION_CREDENTIALS=/path/your/file/creds.json

You are almost set to use the GCP platform. The last remaining thing you need to do is set the configuration file in `conf/gcp.conf` and set your project information, the location and type of machines you want to use:

- `project`: name of your project on GCP
- `zone`: Location of the cloud instances. See `GCP documentation <https://cloud.google.com/compute/docs/regions-zones/?hl=en>`_ for a list of available locations.
- `machineType`: specs of the instance you wish to use. A list of available instances is available on the `GCP documentation <https://cloud.google.com/compute/all-pricing>`_.

For more details about the GCP configuration, see the `nextflow documentation <https://www.nextflow.io/docs/latest/google.html>`_

Configure for Amazon Web Service
--------------------------------

Not implemented yet. If you are interested to set it, have a look at the `nextflow documentation <https://www.nextflow.io/docs/latest/awscloud.html>`_
