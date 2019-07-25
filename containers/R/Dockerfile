FROM r-base:3.6.1
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y nano emacs wget curl

# Install R libraries
RUN Rscript -e "install.packages(c('BiocManager','stringr','seqinr','ggplot2','remotes'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install('dada2')"
RUN Rscript -e "remotes::install_github('tobiasgf/lulu')"

COPY . .

CMD /bin/bash
