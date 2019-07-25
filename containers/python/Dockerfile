FROM python:3.7
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

# Install python libraries
RUN pip3 install ipython biopython numpy pandas itsxpress

# HMMER
RUN curl -L http://eddylab.org/software/hmmer/hmmer-3.2.tar.gz \
	| tar xz
RUN cd hmmer-3.2 && ./configure && make && make install && cd .. && rm -rf hmmer*

# BBtools
RUN curl -L https://sourceforge.net/projects/bbmap/files/latest/download \
	| tar xz -C /usr/local/bin

# VSEARCH
RUN curl -L "https://github.com/torognes/vsearch/releases/download/v2.11.1/vsearch-2.11.1-linux-x86_64.tar.gz" \
	| tar xz -C /usr/local/bin

# Set PATH
ENV PATH /usr/local/bin/vsearch-2.11.1-linux-x86_64/bin:/usr/local/bin/bbmap:PATH

COPY . .

CMD /bin/bash
