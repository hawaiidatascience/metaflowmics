FROM r-base:3.6.1
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y nano wget curl procps
RUN apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev libssh2-1-dev libcairo2-dev libgeos-dev libudunits2-0 libudunits2-dev libgdal-dev

# Install R libraries
RUN Rscript -e "install.packages(c('BiocManager','stringr','seqinr','ggplot2','remotes','doParallel'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install('phyloseq')"
RUN Rscript -e "remotes::install_github('benjjneb/dada2', ref='v1.14')"
RUN Rscript -e "remotes::install_github('tobiasgf/lulu')"
RUN Rscript -e "install.packages(c('stringi', 'igraph', 'vegan'))"

WORKDIR /workspace

COPY . /workspace

CMD /bin/bash
