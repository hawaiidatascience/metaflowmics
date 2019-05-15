### Install conda

Install conda as shown at [bioconda.github.io](https://bioconda.github.io/) :

1. Download [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (python3 version preferred).
2. Install it by running the downloaded script (e.g. `bash Miniconda3-latest-Linux-x86_64.sh`).
3. Select a directory for the installation (e.g. `/opt/conda`)
4. Answer `yes` to the question:
```
Do you wish the installer to initialize Miniconda3
by running conda init? [yes|no]
```
5. Exit the current session and login again to activate conda
6. Modify the channel list as follows:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels r
conda config --add channels conda-forge
```

### Download the workflows

Download the workflows from github:
```
git clone git@github.com:hawaiidatascience/nextflow_cmaiki.git
```

### Create the conda environment "cmaiki"

All the dependencies of the workflow will be installed in the conda environment `cmaiki` with:
```
source build_conda_env`
```

### Run the test for 16S

Run a test with:
```
make test16S-conda 
```
