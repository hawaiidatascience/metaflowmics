FROM ubuntu:18.04
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y wget unzip curl nano procps libreadline7

# Install mothur version
RUN wget "https://github.com/mothur/mothur/releases/download/v.1.43.0/Mothur.Ubuntu_18.zip"
RUN unzip Mothur.Ubuntu_18.zip \
	&& rm -rf __MACOSX Mothur.Ubuntu_18.zip\
	&& mv mothur /usr/local/bin

# Setup PATH
ENV PATH /usr/local/bin/mothur:$PATH

WORKDIR /workspace
COPY . /workspace

CMD /bin/bash
