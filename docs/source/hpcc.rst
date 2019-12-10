.. _hpc_conf:

Configure for a HPCC
--------------------

HPCC require a few more specific parameters to be set:

- `executor`: The HPCC scheduler. Default is SLURM, but many others are available. For a comprehensive list see `nextflow documentation <https://www.nextflow.io/docs/latest/executor.html>`_
- `queue`: List of queues where you wish to submit the processes. Multiple queues are separated with commas.
- `module`: Singularity module location.

