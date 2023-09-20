# Base image for superFreq/v1.4.5,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04

# Depedencies of superFreq:
#  - https://github.com/ChristofferFlensburg/superFreq#dependencies
# Overview:
#  -+ samtools/1.X              # apt-get, installs: 1.13-4
#  -+ bedtools                  # apt-get, installs: 2.30.0
#  -+ vcftools                  # apt-get, installs: 0.1.16-3
#  -+ bcftools                  # apt-get, installs: 1.13-1
#  -+ python3                   # apt-get, installs: 3.10.6-1
#  -+ R>=4.X                    # apt-get, installs: 4.1.2-1
#     @requires: WriteXLS (>= 4.0), BiocGenerics (>= 0.26), 
#                GenomicRanges (>= 1.32.7), IRanges (>= 2.14.12), 
#                S4Vectors (>= 0.18.3), biomaRt (>= 2.36.1), 
#                Rsamtools (>= 1.32.3), R.oo (>= 1.22), 
#                Rsubread (>= 1.30.9), limma (>= 3.36.5), 
#                parallel, MutationalPatterns (>= 1.6.1), 
#                BSgenome.Mmusculus.UCSC.mm10 (>= 1.4.0), 
#                BSgenome.Hsapiens.UCSC.hg38 (>= 1.4.1)
LABEL maintainer="OpenOmics" \
    base_image="ubuntu:22.04" \
    version="v0.1.0"   \
    software="superFreq/v1.4.5" \
    about.summary="superFreq: An R package that analyses cancer exomes" \
    about.home="https://github.com/ChristofferFlensburg/superFreq" \
    about.documentation="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007603" \
    about.tags="Cancer Genomics"

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
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
        # bcftools/1.13-1
        bcftools \ 
        # bedtools/2.30.0
        bedtools \
        build-essential \
        cmake \
        cpanminus \
        curl \
        gawk \
        git \
        gzip \
        libcurl4-openssl-dev \
        libbz2-dev \
        liblzma-dev \
        libmysqlclient-dev \
        libssl-dev \
        libxml2-dev \
        locales \
        # perl/5.34.0-3
        perl \
        pkg-config \
        # python/3.10.6
        python3 \
        python3-pip \
        # R/4.1.2-1
        r-base \
        # devtools/2.4.3
        r-cran-devtools \
        # bioc/1.30.16
        r-cran-biocmanager \
        r-cran-optparse \
        r-cran-tidyverse \
        # genomeinfodb/1.30.1
        r-bioc-genomeinfodb \
        # genomicfeatures/1.46.4
        r-bioc-genomicfeatures \
        # variantannotation/1.40.0
        r-bioc-variantannotation \
        # samtools/1.13
        samtools \
        # vcftools/0.1.16
        vcftools \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
# Perl fix issue
RUN cpanm FindBin Term::ReadLine


############### INSTALL #################
# Install bioinformatics tools/packages,
# Installing useful python packages 
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN pip3 install --upgrade pip \
    && pip3 install matplotlib \
    && pip3 install numpy \
    && pip3 install pandas \
    && pip3 install pybedtools \
    && pip3 install pysam \
    && pip3 install scipy


############ superFreq/1.4.5 ############
# Install superFreq/1.4.5, for more info:
# github.com/ChristofferFlensburg/superFreq/
# setting the timeout to 15 mins due to very
# large size (>600MB) of BSgenome genomes.
RUN Rscript -e \
    'options(timeout = 900); devtools::install_github("ChristofferFlensburg/superFreq", Ncpus = 4)' \
    # Check if installed properly,
    # sometimes the BSgenome package
    # downloads will time out.
    && Rscript -e 'library("superFreq")'

################ VEP/110 ################
# superFreq can also optionally run VEP
# to annotate its results. We are making
# use of this feature. It is important
# to note that the VEP cache dirs are
# not installed within this docker image
# due to size constraints; however, in
# VEP/109 a new option was added to point
# to the plugin path and cache dir path
# with environment variable. This is very
# import because superFreq does not have 
# a built-in method to set/specify this
# location at run time. As so, before 
# running superFreq these paths NEED to
# be export aprori to run time. For more
# information please see the following: 
# https://github.com/Ensembl/ensembl-vep/pull/1267
RUN mkdir -p /opt2/VEP/110/ /opt2/VEP/110/Plugins \
    && wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/110.1.tar.gz -O /opt2/VEP/110/110.1.tar.gz \
    && tar -zvxf /opt2/VEP/110/110.1.tar.gz -C /opt2/VEP/110/ \
    && rm -f /opt2/VEP/110/110.1.tar.gz \
    && cd /opt2/VEP/110/ensembl-vep-release-110.1 \
    # Install Perl dependencies:
    # http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#requirements
    && cpanm DBI DBD::mysql Archive::Zip \
    # Installs API + Faidx/htslib
    && ./INSTALL.pl --NO_UPDATE --AUTO a

ENV PATH="${PATH}:/opt2/VEP/110/ensembl-vep-release-110.1"
WORKDIR /opt2 


################ POST #################
# Add Dockerfile and export environment 
# variables and update permissions
ADD Dockerfile /opt2/superfreq_1-4-5.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2:$PATH"
WORKDIR /data2