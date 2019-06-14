FROM ubuntu:18.04

WORKDIR /code

RUN apt-get update

RUN apt-get install -y wget openjdk-8-jre python3-pip

ENV DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/Pacific/Honolulu /etc/localtime

# Install R
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update && apt-get install -y r-base

RUN wget -qO- https://get.nextflow.io | bash

RUN pip3 install -r ipython biopython pandas fastqc itsxpress

# Install mothur
RUN wget "https://github.com/mothur/mothur/releases/download/v1.42.1/Mothur.linux_64.zip"
RUN unzip "Mothur.linux_64.zip && rm Mothur.linux_64.zip" && rm -rf "__MACOSX"

RUN wget "https://github.com/torognes/vsearch/releases/download/v2.11.1/vsearch-2.11.1-linux-x86_64.tar.gz"
RUN tar -xvzf "vsearch-2.11.1-linux-x86_64.tar.gz"

# Install R libraries
RUN Rscript -e "install.packages(c('stringr','ggplot2','seqinr'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite('dada2')"
RUN Rscript -e "devtools::install_github('tobiasgf/lulu')"

COPY Pipeline-ITS /code/Pipeline-ITS

ENV PATH /mothur:/vsearch-2.11.1-linux-x86_64/bin/:$PATH

CMD /bin/bash
