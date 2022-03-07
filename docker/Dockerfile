FROM rocker/r-ver:3.5.2

LABEL maintainer "Diego Garrido-Martín <diego.garrido@crg.eu>" \
      version "1.0.0" \
      description "Docker image for mvgwas-nf"

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# R:devtools dependencies
       libcurl4-gnutls-dev \
       libssl1.0-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
       libgit2-dev \
       pkg-config \
# R:car
       libnlopt-dev \
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
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.2.3.tar.gz", repos=NULL, type="source");'\
'install.packages(c("devtools", "optparse", "data.table", "seqminer", "vegan", "car", "R.utils"), repos="http://cloud.r-project.org/")'

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

## Install MANTA
RUN R -e 'devtools:::install_github("dgarrimar/manta")'

## Change ENTRYPOINT
ENTRYPOINT /bin/bash

