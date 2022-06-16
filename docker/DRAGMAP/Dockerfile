# Dragmap Dockerfile
# Using Ubuntu Focal (20.04 LTS) as base image
FROM ubuntu:20.04
LABEL maintainer=kuhnsa@nih.gov

# Overview of Dependencies
#  • C++17 compatible compiler (e.g gcc-c++ >= 7.1)
#  • GNU make >= 3.82 and cmake
#  • Compression Library: zlib
#  • Boost library: libboost-all-dev
#  • For unit tests: googletest (>= v1.6) and libgtest-dev
#  • Git: get version information

# Define release with --build-arg VERSION=$RELEASE_VERSION
# where the default = 1.3.0
ARG VERSION=1.3.0

# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with host filesyetem 
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        libboost-all-dev \
        libgtest-dev \
        libz-dev \
        g++ \
        gcc \
        git \
        googletest \
        locales \
        cmake \
        make \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Install DRAGMAP from official release
RUN wget --no-check-certificate https://github.com/Illumina/DRAGMAP/archive/refs/tags/${VERSION}.tar.gz \
    && tar -xvf ${VERSION}.tar.gz \
    && rm ${VERSION}.tar.gz \
    && cd DRAGMAP-${VERSION} \
    && make \
    && make install \
    && dragen-os --help

WORKDIR /opt2 

# Add Dockerfile and export env variables
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2