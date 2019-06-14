FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/Pacific/Honolulu /etc/localtime

RUN apt-get update 
RUN apt-get install -y wget openjdk-8-jre python3-pip libcurl4-openssl-dev libxml2-dev libssl-dev

# Install R
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update && apt-get install -y r-base

RUN wget -qO- https://get.nextflow.io | bash

# Install R libraries
RUN Rscript -e "install.packages(c('stringr','ggplot2','seqinr','BiocManager'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install('dada2')"
RUN Rscript -e "devtools::install_github('tobiasgf/lulu')"

RUN apt-get install -y nano emacs25 vim

# Install vll
RUN wget -qO- "https://get.haskellstack.org" | sh
RUN wget "https://github.com/w9/vll-haskell/releases/download/v0.0.1-hotfix2/vl" && chmod +x vl
RUN wget "https://github.com/w9/vll-haskell/releases/download/v0.0.1-hotfix2/vll" && chmod +x vll

RUN pip3 install ipython biopython pandas itsxpress

# Install VSEARCH
RUN wget "https://github.com/torognes/vsearch/releases/download/v2.11.1/vsearch-2.11.1-linux-x86_64.tar.gz"
RUN tar -xvzf "vsearch-2.11.1-linux-x86_64.tar.gz" && rm -f "vsearch-2.11.1-linux-x86_64.tar.gz"

# Install mothur
RUN wget "https://github.com/mothur/mothur/releases/download/v1.42.1/Mothur.linux_64.zip"
RUN unzip "Mothur.linux_64.zip" && rm -f "Mothur.linux_64.zip" && rm -rf "__MACOSX"

# Copy databases
WORKDIR /code
COPY databases.tar.gz /code/databases.tar.gz
RUN tar -xvzf /code/databases.tar.gz && rm -rf /code/databases.tar.gz

# Setup PATH
RUN mv /mothur /vl /vll "/vsearch-2.11.1-linux-x86_64" /usr/local/bin/
ENV PATH /usr/local/bin/mothur:/usr/local/bin/vsearch-2.11.1-linux-x86_64/bin:$PATH

# Copy files
COPY Pipeline-ITS /code/Pipeline-ITS
COPY Pipeline-16S /code/Pipeline-16S

CMD /bin/bash
