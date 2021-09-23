from debian:bookworm

MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

RUN apt-get update -y && apt-get install -y wget
RUN wget https://github.com/rcedgar/muscle/raw/main/binaries/muscle5.0.1278_linux64 \
    && mv muscle5.0.1278_linux64 /usr/local/bin/muscle \
    && chmod +x /usr/local/bin/muscle

WORKDIR /workspace

CMD /bin/bash
