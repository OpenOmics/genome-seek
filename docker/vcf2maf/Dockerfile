# vcf2maf Dockerfile
# Using VEP 106 from Ensembl as base image, 
# which uses @ubuntu/18.04 as its base image: 
# hub.docker.com/r/ensemblorg/ensembl-vep/tags
# VEP/106 matches GENCODE/v40 for GRCh38.
# @Dockerfile of base image: 
# https://github.com/Ensembl/ensembl-vep/blob/release/106/docker/Dockerfile
FROM ensemblorg/ensembl-vep:release_106.0

LABEL maintainer=kuhnsa@nih.gov

# Overview of Dependencies
# https://github.com/mskcc/vcf2maf/tree/v1.6.21

# Define release with --build-arg VERSION=$RELEASE_VERSION
# where the default = 1.6.21
ARG VERSION=1.6.21

# Create Container filesystem specific 
# working directory and opt directories
USER root
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system 
# packages required for your 
# project. If you need extra 
# system packages add them here.
# argparse.bash requires argparse
# python package.
# Installs python/3.6, which
# includes argparse by default
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        locales \
        libncurses5-dev \
        libncursesw5-dev \
        libpng-dev \
        libssl-dev \
        libcurl4-gnutls-dev \
        libmysqlcppconn-dev \
        liblzma-dev \
        libbz2-dev \
        openssl \
        python3 \
        python3-pip \
        rsync \
        uuid-dev \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Install htslib/1.10
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 \
        && tar -vxjf htslib-1.10.tar.bz2 \
        && rm htslib-1.10.tar.bz2 \
        && cd htslib-1.10 \
        && make \
        && cd /opt2

# Install samtools/1.10
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
        && tar -vxjf samtools-1.10.tar.bz2 \
        && rm samtools-1.10.tar.bz2 \
        && cd samtools-1.10 \
        && make \
        && cd /opt2

# Install bcftools/1.10
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.10/bcftools-1.10.tar.bz2 \
        && tar -vxjf bcftools-1.10.tar.bz2 \
        && rm bcftools-1.10.tar.bz2 \
        && cd bcftools-1.10 \
        && make \
        && cd /opt2

# Install UCSC/377
# vcf2maf needs UCSC liftOver
RUN wget --no-check-certificate http://hgdownload.cse.ucsc.edu/admin/exe/userApps.archive/userApps.v377.src.tgz \
        && tar -xvzf userApps.v377.src.tgz \
        && rm userApps.v377.src.tgz \
        && cd userApps \
        && make \
        && find /opt2/userApps/bin/ -type f -not -iname 'lift*' -delete \
        && rm -rf /opt2/userApps/kent/src/ \
        && cd /opt2

# Install vcf2maf/v1.6.21
RUN wget --no-check-certificate https://github.com/mskcc/vcf2maf/archive/refs/tags/v${VERSION}.tar.gz \
        && tar -xvf v${VERSION}.tar.gz \
        && rm v${VERSION}.tar.gz \
        && chmod a+rx /opt2/vcf2maf-${VERSION}/*.pl

# Add build Dockerfile to /opt2
# and export environment variables
# and ensure proper permissions
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 \
    && ln -sf /usr/bin/python3 /usr/bin/python
ENV PATH="$PATH:/opt2/bcftools-1.10:/opt2/samtools-1.10:/opt2/htslib-1.10:/opt2/userApps/bin:/opt2:/opt2/vcf2maf-${VERSION}"
ENV VEP_HOME="/opt/vep/src/ensembl-vep"

# Reset working directory
WORKDIR /data2