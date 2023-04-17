# Base image for mr-seek,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="genome-seek" \
    about.summary="An awesome Clinical Whole Genome Sequencing Pipeline" \
    about.home="https://github.com/OpenOmics/genome-seek" \
    about.documentation="https://openomics.github.io/genome-seek/" \
    about.tags="Variant-calling"


############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with host filesyetem 
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 
# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone


############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apt-transport-https \
        autotools-dev \
        automake \
        autoconf \
        bc \
        bcftools \
        build-essential \
        bzip2 \
        ca-certificates \
        cmake \
        cpanminus \
        curl \
        dirmngr \
        g++ \
        gawk \
        gcc \ 
        gfortran \
        git \
        gnupg \
        gzip \
        libcairo2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev  \
        libgmp-dev \
        libtool \
        libz-dev \
        locales \
        make \
        openjdk-8-jdk \
        pandoc \
        parallel \
        pigz \
        python2 \
        python3 \
        python3-pip \
        software-properties-common \
        unzip \
        wget \
        zlib1g \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8


############### INSTALL ################
# Make Python3 the default interpreter
# and install python packages
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
	  && pip3 install argparse \
	  && pip3 install numpy \
      && pip3 install pysam \
	  && pip3 install scipy

# Perl fix issue
RUN cpanm FindBin Term::ReadLine

# Install R (4.2.2) and R packages.
# For more information, check out: 
# https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-get update -qq \
    # Add the signing key r-project repos
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install \
        libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
# Install CRAN packages
RUN Rscript -e 'install.packages(c("data.table", "dplyr", "optparse", "remotes", "devtools"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("R.utils", "tidyverse", "ggpubr", "ggsignif", "rstatix")), repos="http://cran.r-project.org")'
# Install Bioconductor packages 
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("limma", "edgeR", "rtracklayer", "maftools"))'

# TODO: apt-get common bioinformatics tools
# ...

# TODO: manually install/build tools from Github 
# ...


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek.dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2