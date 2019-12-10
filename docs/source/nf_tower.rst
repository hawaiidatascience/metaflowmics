.. _nf_tower:

Monitoring your runs with Nextflow Tower
========================================

`Nextflow Tower <https://tower.nf/>`_ is a monitoring tool for your nextflow runs. Through a web interface, it tracks all processes launched by your nextflow run and returns some useful statistics on your jobs, such as the CPU and RAM usage or the running time.

To include this monitoring tool in your runs, you will need to:

- Create sign in `nextflow tower website <https://tower.nf/>`_
- Generate a token:

  #. Once you are logged in, click on `Your tokens` on the top left of the screen
  #. Click on `New token`

- Provide this token to nextflow (as described in `tower documentation <https://tower.nf/welcome>`_) by either:
  - Setting the environment variable `TOWER_ACCESS_TOKEN` to your token value
  - Setting the nextflow configuration variable, `tower.accessToken` to your token value

Then, you can simply run your nextflow pipelines with the `with-tower` flag. The summaries will appear on your personal account on `Nextflow Tower website <https://tower.nf>`_
