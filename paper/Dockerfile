FROM ubuntu:20.04

ARG varnote_version=1.2.0
ARG echtvar_version=0.1.3
ARG slivar_version=0.2.7
ARG htslib_version=1.14

ENV DEBIAN_FRONTEND noninteractive
ENV LC_ALL C

RUN \
    # install packages dependencies
    apt-get update -yqq && \
    apt-get install -yqq \
        curl time \
        r-base r-base-dev \
        git \
        bash \
        openjdk-17-jre-headless \
        locales \
        wget \
        libbz2-dev zlib1g-dev \
        liblzma-dev \
    && apt-get clean && \
    \
    # configure locale, see https://github.com/rocker-org/rocker/issues/19
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8

RUN wget -q https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 && rm *.tar.bz2 && cd htslib* && \
    ./configure && make -j4 install && \
    cd .. && \
    wget -q https://github.com/samtools/bcftools/releases/download/${htslib_version}/bcftools-${htslib_version}.tar.bz2 \
    && tar xjf bcftools-${htslib_version}.tar.bz2 && rm *.tar.bz2 && cd bcftools* && \
    ./configure && make -j4 install

RUN mkdir -p /opt/programs/ && \
    wget -q https://github.com/mulinlab/VarNote/releases/download/v${varnote_version}/VarNote-${varnote_version}.zip && \
    unzip VarNote-${varnote_version}.zip && mv VarNote-${varnote_version}.jar /opt/programs/VarNote.jar && \
    rm -f VarNote* && \
    wget -qO /usr/local/bin/echtvar https://github.com/brentp/echtvar/releases/download/v${echtvar_version}/echtvar && \
    chmod +x /usr/local/bin/echtvar && \
    wget -qO /usr/local/bin/slivar https://github.com/brentp/slivar/releases/download/v${slivar_version}/slivar && \
    chmod +x /usr/local/bin/slivar && \
    echo '#!/bin/sh\njava -Xmx5G -jar /opt/programs/VarNote.jar "$@"' > /usr/local/bin/varnote && \
    chmod +x /usr/local/bin/varnote
