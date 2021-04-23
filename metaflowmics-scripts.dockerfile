FROM continuumio/miniconda:latest
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

RUN conda install -y -c conda-forge -c bioconda -c r \
    r-stringr r-seqinr bioconductor-dada2 r-remotes \
    r-dplyr r-tidyr r-data.table

RUN Rscript -e "remotes::install_github('tobiasgf/lulu')"

WORKDIR /workspace

CMD /bin/bash
