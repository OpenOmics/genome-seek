
# Base image for genome-seek,
# uses Ubuntu Jammy (LTS)
FROM ubuntu:22.04


############# OVERVIEW #############
#  • python+pip3/3.10 (apt-get python3)
#    - argparse  (standard library >= 3.5)
#    - pandas
#    - numpy
#    - scipy
#  • R/4.2.2
#    - maftools
#  • bc, bzip2, curl, gawk, git, gzip, 
#      java8, parallel, pandoc, pigz, 
#      python2.7, unzip, wget
#  • gatk3, gatk4, picard
#  • bcftools, vcftools, samtools
#  • deepvariant, muse, strelka
#  • fastp, bwa_mem2, samblaster
#  • circos, canvas, peddy, manta
#  • hla_la
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
        libcurl4-gnutls-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev  \
        libgmp-dev \
        libgsl0-dev \
        liblzma-dev \
        libncurses5-dev \
        libperl-dev \
        libssl-dev \
        libtool \
        libz-dev \
        libbz2-dev \
        locales \
        make \
        openjdk-8-jdk \
        pandoc \
        parallel \
        perl \
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

# Install R (4.2.3) and R packages.
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
RUN Rscript -e 'install.packages(c("data.table", "dplyr", "optparse", "remotes", "devtools", "plotly"), Ncpus = 4, repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("R.utils", "tidyverse", "ggpubr", "ggsignif", "rstatix"), Ncpus = 4, repos="http://cran.r-project.org")'
# Install Bioconductor packages 
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("limma", "edgeR", "rtracklayer", "maftools"), Ncpus = 4)'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("copynumber", "VariantAnnotation"), Ncpus = 4)'

# Install common bioinformatic 
# tools with apt-get
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        # bedtools 2.30.0+dfsg-2
        bedtools \
        # bowtie 1.3.1-1
        bowtie \
        # bowtie2 2.4.4-1
        bowtie2 \
        # samblaster 0.1.26-3
        samblaster \
        # tabix 1.13+ds-2build1
        tabix \
        # circos 0.69.9+dfsg-2
        circos \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Manually install/build tools from source
# Install GATK4 (GATK/4.3.0.0)
# Requires Java8 or 1.8
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip \
    && unzip /opt2/gatk-4.3.0.0.zip \
    && rm /opt2/gatk-4.3.0.0.zip \
    && /opt2/gatk-4.3.0.0/gatk --list
ENV PATH="/opt2/gatk-4.3.0.0:$PATH"

# Install last release of GATK3 (GATK/3.8-1)
# Only being used for the CombineVariants
# command that is not available in GATK4
# Available via env variable: $GATK_JAR
# Requires Java8 or 1.8
RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    && tar -xvjf /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    && rm /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
ENV GATK_JAR="/opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

# Install picard/2.27.5,
# requires java8 or 1.8
RUN mkdir -p /opt2/picard/2.27.5/ \
    && wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -P /opt2/picard/2.27.5/
ENV PICARDJARPATH=/opt2/picard/2.27.5/ \
    PICARD_JARPATH=/opt2/picard/2.27.5/ \
    PICARDJAR=/opt2/picard/2.27.5/picard.jar \
    PICARD_JAR=/opt2/picard/2.27.5/picard.jar \
    PICARD=/opt2/picard/2.27.5/picard.jar

# Install samtools/1.17 + htslib,
# see the following for more info:
# https://github.com/samtools/samtools/blob/develop/INSTALL
RUN mkdir -p /opt2/samtools/1.17/ \
    && wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 -P /opt2/samtools/1.17/ \
    && tar -xjvf /opt2/samtools/1.17/samtools-1.17.tar.bz2 -C /opt2/samtools/1.17/\
    && rm /opt2/samtools/1.17/samtools-1.17.tar.bz2 \
    && ls -larth /opt2/samtools/1.17/ \
    && cd /opt2/samtools/1.17/samtools-1.17/ \
    && ./configure \
    && make \
    && make install
WORKDIR /opt2

# Install bcftools/1.17 + htslib,
# see the following for more info:
# https://github.com/samtools/bcftools/blob/develop/INSTALL
RUN mkdir -p /opt2/bcftools/1.17/ \
    && wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 -P /opt2/bcftools/1.17/ \
    && tar -xjvf /opt2/bcftools/1.17/bcftools-1.17.tar.bz2 -C /opt2/bcftools/1.17/ \
    && rm /opt2/bcftools/1.17/bcftools-1.17.tar.bz2 \
    && cd /opt2/bcftools/1.17/bcftools-1.17/ \
    && ./configure \
    && make \
    && make install
WORKDIR /opt2

# Install vcftools/0.1.16,
# see the following for more info:
# https://vcftools.github.io/examples.html
RUN mkdir -p /opt2/vcftools/0.1.16/ \
    && wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz -P /opt2/vcftools/0.1.16/ \
    && tar -vxzf /opt2/vcftools/0.1.16/vcftools-0.1.16.tar.gz -C /opt2/vcftools/0.1.16/ \
    && rm /opt2/vcftools/0.1.16/vcftools-0.1.16.tar.gz \
    && export PERL5LIB=/opt2/vcftools/0.1.16/vcftools-0.1.16/src/perl/ \
    && cd /opt2/vcftools/0.1.16/vcftools-0.1.16/ \
    && ./configure \
    && make \
    && make install
WORKDIR /opt2


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek.dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2