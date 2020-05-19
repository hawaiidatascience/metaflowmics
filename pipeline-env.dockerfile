FROM ubuntu:18.04
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y apt-transport-https software-properties-common wget unzip nano emacs procps default-jre

#-------------------------------------------------#
#               Mothur requirements               #
#-------------------------------------------------#
RUN apt-get install -y libreadline7

#-------------------------------------------------#
#              R package requirements             #
#-------------------------------------------------#

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt-get update && apt-get install -y r-base libcurl4-openssl-dev libxml2-dev libssl-dev libssh2-1-dev libcairo2-dev libgeos-dev libudunits2-0 libudunits2-dev libgdal-dev

#-------------------------------------------------#
#            Requirements for ITSxpress           #
#       (VSEARCH already provided by Mothur)      #
#-------------------------------------------------#

# HMMER
RUN wget -qO- http://eddylab.org/software/hmmer/hmmer-3.2.tar.gz \
	| tar xz
RUN cd hmmer-3.2 && ./configure && make && make install && cd .. && rm -rf hmmer*

# BBtools
RUN wget -qO- https://sourceforge.net/projects/bbmap/files/latest/download \
	| tar xz -C /usr/local/bin

#-------------------------------------------------#
#                   R libraries                   #
#-------------------------------------------------#

RUN Rscript -e "install.packages(c('BiocManager','stringr','seqinr','ggplot2','remotes','doParallel', 'stringi', 'igraph', 'vegan'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install('phyloseq')"
RUN Rscript -e "remotes::install_github('benjjneb/dada2', ref='v1.16')"
RUN Rscript -e "remotes::install_github('tobiasgf/lulu')"

#-------------------------------------------------#
#                  python packages                #
#-------------------------------------------------#

RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install ipython biopython numpy pandas itsxpress matplotlib seaborn h5py scikit-learn bokeh scipy

#-------------------------------------------------#
#                   Mothur v1.43                  #
#-------------------------------------------------#

RUN wget "https://github.com/mothur/mothur/releases/download/v.1.44.1/Mothur.Ubuntu_18.zip"
RUN unzip Mothur.Ubuntu_18.zip \
	&& rm -rf __MACOSX Mothur.Ubuntu_18.zip\
	&& mv mothur /usr/local/bin

#-------------------------------------------------#
#                      FastQC                     #
#-------------------------------------------------#

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
	&& unzip fastqc_v0.11.8.zip && rm fastqc_v0.11.8.zip \
    && chmod a+x FastQC/fastqc \
	&& mv FastQC /usr/local/bin

#-------------------------------------------------#
#    Environment variables and work directoy      #
#-------------------------------------------------#

ENV PATH /usr/local/bin/mothur:/usr/local/bin/bbmap:/usr/local/bin/FastQC:$PATH

WORKDIR /workspace
COPY . /workspace

CMD /bin/bash
