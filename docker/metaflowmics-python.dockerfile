FROM continuumio/miniconda:latest
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

RUN apt-get update --allow-releaseinfo-change -y && apt-get install -y procps

RUN conda install -y -c conda-forge -c bioconda \
    python">=3.6" pandas">=1" holoviews=1.14.5 datatable scipy biopython h5py scikit-learn

WORKDIR /workspace

CMD /bin/bash
