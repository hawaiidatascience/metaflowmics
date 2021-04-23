FROM continuumio/miniconda:latest
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

# RUN apt-get update && apt-get install -y apt-transport-https software-properties-common wget unzip nano emacs procps default-jre gcc

#-------------------------------------------------#
#            Requirements for ITSxpress           #
#       (VSEARCH already provided by Mothur)      #
#-------------------------------------------------#

RUN conda install -y -c conda-forge -c bioconda -c r \
    r-stringr r-stringi r-seqinr r-ggplot2 r-doparallel r-igraph r-vegan r-remotes \
    bioconductor-phyloseq "bioconductor-dada2 >=1.18" \
    itsxpress "fastqc>=0.11.9" "python>=3.6" pip

RUN conda install -y -c bioconda -c conda-forge r-remotes r-dplyr vsearch
RUN Rscript -e "remotes::install_github('tobiasgf/lulu')"

#-------------------------------------------------#
#                  python packages                #
#-------------------------------------------------#

RUN pip3 install ipython biopython numpy pandas itsxpress matplotlib seaborn h5py scikit-learn bokeh scipy

#-------------------------------------------------#
#    Environment variables and work directory      #
#-------------------------------------------------#

WORKDIR /workspace
COPY . /workspace

CMD /bin/bash
