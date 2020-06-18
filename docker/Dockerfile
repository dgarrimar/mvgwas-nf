FROM rocker/r-ver:3.5.2

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# devtools dependencies
       libcurl4-openssl-dev \
       libssl-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
# bcftools dependencies
       libbz2-dev \
       liblzma-dev \
# nextflow trace (ps)
       procps \
# other
       wget \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install R packages
RUN R -e 'install.packages(c("devtools", "optparse", "data.table", "seqminer", "vegan"), repos="http://cloud.r-project.org/")'

## Install mlm
RUN R -e 'devtools:::install_github("dgarrimar/mlm")'

## Change ENTRYPOINT
ENTRYPOINT /bin/bash

## Install bcftools and HTSlib
RUN wget --quiet https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
    && tar jxf bcftools-1.10.2.tar.bz2 \
    && cd bcftools-1.10.2 \
    && ./configure \
    && make \
    && make install

RUN cd bcftools-1.10.2/htslib-1.10.2 \
    && ./configure \
    && make \
    && make install
