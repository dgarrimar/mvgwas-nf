FROM rocker/r-ver:3.5.2

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# devtools dependencies
       libcurl4-openssl-dev \
       libssl1.0-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
# libgit2 
       cmake \
       pkg-config \
       git \
# car
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

## Install libgit2 (compile from source as in https://github.com/libgit2/pygit2/issues/999#issuecomment-614510556) 
RUN git clone --depth=1 -b v1.0.0 https://github.com/libgit2/libgit2.git ~/libgit2_src \
    && cd ~/libgit2_src \
    && cmake . -DBUILD_CLAR=OFF -DCMAKE_BUILD_TYPE=Release -DEMBED_SSH_PATH=~/libssh2_src -DCMAKE_INSTALL_PREFIX=~/libgit2 \
    && cmake --build . --target install \
    && cp -r ~/libgit2/* /usr/bin \
    && cp -r ~/libgit2/* /usr/local

## Install R packages
RUN R -e 'install.packages("lme4", repos="http://cloud.r-project.org/");'\
'install.packages("https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-7.tar.gz", repos=NULL, type="source");'\
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

