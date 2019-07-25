FROM ubuntu:18.04
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y nano emacs wget curl python3-pip
RUN apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev libssh2-1-dev libcairo2-dev libgeos-dev libudunits2-0 libudunits2-dev libgdal-dev locales

RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# Install libreadline6 for Mothur
RUN echo "deb http://archive.ubuntu.com/ubuntu/ xenial main" >> /etc/apt/sources.list \
	&& apt-get update \
	&& apt-get install -y libreadline6

# Install R
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update && apt-get install -y r-base

# Install R libraries
RUN Rscript -e "install.packages(c('BiocManager','stringr','seqinr','ggplot2'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install('dada2')"
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "remotes::install_github('tobiasgf/lulu')"

# Install vll
RUN wget -qO- "https://get.haskellstack.org" | sh
RUN wget -o /usr/local/bin/vl "https://github.com/w9/vll-haskell/releases/download/v0.0.1-hotfix2/vl"
RUN wget -o /usr/local/bin/vll "https://github.com/w9/vll-haskell/releases/download/v0.0.1-hotfix2/vll"
RUN chmod +x /usr/local/bin/vl*

RUN useradd -ms /bin/bash developer
WORKDIR /home/developer

# Install BBtools
RUN curl -L https://sourceforge.net/projects/bbmap/files/latest/download \
	| tar xz -C /usr/local/bin

# Install hmmer
RUN curl -L http://eddylab.org/software/hmmer/hmmer-3.2.tar.gz \
	| tar xz
RUN cd hmmer-3.2 && ./configure && make && make install && cd .. && rm -rf hmmer*

# Install fastx toolkit
RUN apt-get install -y fastx-toolkit

# Download databases
RUN mkdir databases
RUN curl -L "https://mothur.org/w/images/3/32/Silva.nr_v132.tgz" | tar xz -C databases
RUN wget "https://files.plutof.ut.ee/doi/A1/C9/A1C964DFB03C2A1B37FA16784BA739C88F0941AC68560CEA54DD707F1CF00AC4.zip" -O unite_db.zip \
	&& unzip unite_db.zip && rm unite_db.zip
RUN sed -e 's/|SH.*//g' -e 's/;/,/g' -e 's/__/:/g' -e 's/|/;tax=d:Eukaryota,/g' UNITE_public_all_02.02.2019.fasta \
	| iconv -f utf8 -t ascii//TRANSLIT > databases/unite_all_eukaryotes.fasta

# Install VSEARCH
RUN curl -L "https://github.com/torognes/vsearch/releases/download/v2.11.1/vsearch-2.11.1-linux-x86_64.tar.gz" | tar xz -C /usr/local/bin

# Install mothur versions
RUN wget "https://github.com/mothur/mothur/releases/download/v1.42.0/Mothur.linux_64.zip"
RUN unzip Mothur.linux_64.zip \
	&& rm -rf __MACOSX Mothur.linux_64.zip\
	&& mv mothur /usr/local/bin

RUN wget "https://github.com/mothur/mothur/releases/download/v1.41.3/Mothur.linux_64.zip"
RUN unzip Mothur.linux_64.zip \
	&& rm -rf __MACOSX Mothur.linux_64.zip\
	&& mv mothur /usr/local/bin/mothur_1-41-3

RUN pip3 install ipython biopython pandas itsxpress

RUN apt-get install -y default-jre

RUN wget "https://github.com/nextflow-io/nextflow/releases/download/v19.04.1/nextflow-19.04.1-all"\
	&& chmod 777 nextflow-19.04.1-all\
	&& mv nextflow-19.04.1-all /usr/local/bin/nextflow

# Setup PATH
ENV PATH /usr/local/bin/mothur:/usr/local/bin/vsearch-2.11.1-linux-x86_64/bin:/usr/local/bin/bbmap:$PATH

USER developer

COPY . .

CMD /bin/bash
