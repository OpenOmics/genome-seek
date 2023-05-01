############# OVERVIEW #############
# genome-seek cnv docker image,
# uses genome-seek as a base image:
FROM skchronicles/genome-seek:v0.1.0
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="skchronicles/genome-seek:v0.1.0" \
    version="v0.1.0"   \
    software="genome-seekq cnv" \
    about.summary="An awesome Clinical Whole Genome Sequencing Pipeline" \
    about.home="https://github.com/OpenOmics/genome-seek" \
    about.documentation="https://openomics.github.io/genome-seek/" \
    about.tags="Copy number variation calling"

# Base image already contains:
# bcftools, vcftools, circos, gawk, 
# R/4.X, R-ggplot2, R-tidyverse, 
# R-copynumber, R-VariantAnnotation
# Java 8 is installed in the base 
# image; however hmftools requires
# java 11+.
# @INSTALLS:
#  • java 11+       ✓ (hmftool jars included in resource bundle)
#  • peddy          ✓
#  • dotnet         ✓ (ASP.NET Core 2.1 Runtime)
#  • canvas         ✓ (requires dotnet)
#  • annotsv/3.3.5  ✓ (annotations included in resource bundle)
#  • tcl/tcllib     ✓
WORKDIR /opt2
# Install Java 11+ for hmftools,
# installing dotnet runtime dependencies
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        openjdk-11-jdk \
        libicu70 \
        liblttng-ust1 \
        libssl3 \
        libunwind8 \
        tcl \
        tcllib \
        tcsh \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Installing dotnet for Canvas, for more info see:
# https://github.com/Illumina/canvas/issues/99
# This version of dotnet requires libssl/1.1
RUN mkdir -p /opt2/dotnet/2.1.30/ \
    && wget https://download.visualstudio.microsoft.com/download/pr/d6040f80-8343-4771-9c02-dbc9a35ac88a/68e74e6e46cf36fa1a50f68af6831d6d/aspnetcore-runtime-2.1.30-linux-x64.tar.gz -P /opt2/dotnet/2.1.30/ \
    && tar -vxzf /opt2/dotnet/2.1.30/aspnetcore-runtime-2.1.30-linux-x64.tar.gz -C /opt2/dotnet/2.1.30/ \
    && rm /opt2/dotnet/2.1.30/aspnetcore-runtime-2.1.30-linux-x64.tar.gz \
    && wget https://www.openssl.org/source/openssl-1.1.1c.tar.gz -P /opt2/dotnet/lib \
    && tar -xzvf /opt2/dotnet/lib/openssl-1.1.1c.tar.gz -C /opt2/dotnet/lib \
    && rm /opt2/dotnet/lib/openssl-1.1.1c.tar.gz \ 
    && cd /opt2/dotnet/lib/openssl-1.1.1c/ \
    && ./config \
            --prefix=/opt2/dotnet/lib/openssl \
            --openssldir=/opt2/dotnet/lib/openssl \
    && make \
    && make install
ENV PATH="${PATH}:/opt2/dotnet/2.1.30"
ENV LD_LIBRARY_PATH="/opt2/dotnet/lib/openssl"

############## INSTALL #############
# Install Canvas/1.40.0.1613+master, for more info see:
# https://github.com/Illumina/canvas/
# Creating a wrapper to Canvas.dll that
# is similar to Biowulf HPC setup.
RUN mkdir -p /opt2/canvas/1.40.0.1613.master/ \
    && wget https://github.com/Illumina/canvas/releases/download/1.40.0.1613%2Bmaster/Canvas-1.40.0.1613.master_x64.tar.gz -P /opt2/canvas/1.40.0.1613.master/ \
    && tar -vxzf /opt2/canvas/1.40.0.1613.master/Canvas-1.40.0.1613.master_x64.tar.gz -C /opt2/canvas/1.40.0.1613.master/ \
    && rm /opt2/canvas/1.40.0.1613.master/Canvas-1.40.0.1613.master_x64.tar.gz \
    && rmdir /opt2/canvas/1.40.0.1613.master/1.40.0.1613+master/ \
    && mkdir -p /opt2/canvas/1.40.0.1613.master/bin \
    && echo '#!/bin/bash' > /opt2/canvas/1.40.0.1613.master/Canvas.sh \
    && echo 'dotnet /opt2/canvas/1.40.0.1613.master/Canvas-1.40.0.1613+master_x64/Canvas.dll "$@"' >> /opt2/canvas/1.40.0.1613.master/Canvas.sh \
    && chmod a+rx /opt2/canvas/1.40.0.1613.master/Canvas.sh \
    && ln -sf /opt2/canvas/1.40.0.1613.master/Canvas.sh /opt2/canvas/1.40.0.1613.master/bin/Canvas.dll
ENV PATH="${PATH}:/opt2/canvas/1.40.0.1613.master/bin"

# Install peddy/0.4.8, for more info see:
# https://github.com/brentp/peddy/,
# Checkout the latest commit that is associated with
# the 0.4.8 version of peddy (includes nothing super
# important but better than just running git clone)
RUN mkdir -p /opt2/peddy/0.4.8 \
    && git clone https://github.com/brentp/peddy.git /opt2/peddy/0.4.8/ \
    && cd /opt2/peddy/0.4.8 \
    && git checkout 304aabbac88c9f70e48f312d9121cf9b68c195d0 \
    && pip install -r requirements.txt \
    && pip install --editable .

# Install AnnotSV/3.3.5, for more info see:
# https://lbgi.fr/AnnotSV/downloads
# https://github.com/lgmgeo/AnnotSV
# https://lbgi.fr/AnnotSV/Documentation/README.AnnotSV_3.3.pdf
# @Required libraries:
#  • TCL/Tk                    ✓ Already installed, see above
#  • tcllib (tar/csv packages) ✓ Already installed, see above
#  • bedtools (>=2.25)         ✓ Already in base image
#  • vcftools (>=1.9)          ✓ Already in base image
#  • java (>=8)                ✓ Already in base image
#  • python (>=3.88)           ✓ Already in base image
#    • python-pandas           ✓ Already in base image
#    • python-pyfaidx          - Needs to be installed
#    • python-natsort          - Needs to be installed
#  • Annotations (human)       + Included in resource bundle
# AnnotSV has large resource bundles, the human annotation
# is about 6.6 GB when it is extracted. This is too large 
# to include in the docker image. As so, AnnotSV will be 
# installed in the docker image AND the annotation files 
# will be included in the resource bundle of the pipeline, 
# please see this issue on Github for more information:
# https://github.com/lgmgeo/AnnotSV/issues/177 
RUN mkdir -p /opt2/annotsv/3.3.5/ \
    && pip install natsort pyfaidx \
    && wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.3.5.tar.gz -P /opt2/annotsv/3.3.5/ \
    && tar -vxzf /opt2/annotsv/3.3.5/v3.3.5.tar.gz -C /opt2/annotsv/3.3.5/ \
    && rm /opt2/annotsv/3.3.5/v3.3.5.tar.gz \
    && cd /opt2/annotsv/3.3.5/AnnotSV-3.3.5 \
    && make install
ENV ANNOTSV=/opt2/annotsv/3.3.5/AnnotSV-3.3.5 \
    PATH="${PATH}:/opt2/annotsv/3.3.5/AnnotSV-3.3.5/bin"


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek_cnv.dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2