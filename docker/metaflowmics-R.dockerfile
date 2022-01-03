FROM continuumio/miniconda:latest
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

RUN apt-get update --allow-releaseinfo-change -y && apt-get install -y procps

RUN conda install -y -c conda-forge -c bioconda \
    r">=4.0" r-stringr r-ellipsis">=0.3.2" r-seqinr bioconductor-dada2 r-remotes \
    r-dplyr r-tidyr r-data.table bioconductor-phyloseq">=1.34.0"

RUN Rscript -e "remotes::install_github('tobiasgf/lulu')"

WORKDIR /workspace

CMD /bin/bash
