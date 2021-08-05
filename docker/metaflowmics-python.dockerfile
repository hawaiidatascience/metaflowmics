FROM continuumio/miniconda:latest
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

RUN conda install -y -c conda-forge -c bioconda \
    python">=3.6" pandas">=1" holoviews datatable scipy biopython h5py scikit-learn

WORKDIR /workspace

CMD /bin/bash
