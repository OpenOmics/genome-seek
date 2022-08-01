# Octopus Dockerfile
# Using Ubuntu Jammy (22.04 LTS) as base image
# https://luntergroup.github.io/octopus/docs/installation
FROM ubuntu:22.04
LABEL maintainer=kuhnsa@nih.gov

# Overview of Dependencies
#  • C++14 compatible compiler   (e.g gcc-c++ >= 9.3)
#  • git >= 2.25:                          git (2.34)
#  • Boost library >= 1.65: libboost-all-dev (1.74.0)
#  • htslib >= 1.4 AND != 1.12:     libhts-dev (1.13)
#  • GMP >= 5.1.0:                 libgmp-dev (6.2.1)
#  • GNU make and cmake >= 3.9:        cmake (3.22.1)
#  • python version >= 3:            python3 (3.10.4)
#    • Python Package: distro

# Define release with --build-arg VERSION=$RELEASE_VERSION
# where the default = v0.7.4
ARG VERSION=v0.7.4

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
        bcftools \
        curl \
        git \
        locales \
        libboost-all-dev \
        libgmp-dev \
        libz-dev \
        libhts-dev \
        cmake \
        make \
        python3 \
        python3-pip \
        tabix \
        vcftools \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Install Python Packages
RUN pip3 install --upgrade pip \
	  && pip3 install distro

# Install Octopus from tagged version
RUN git clone --branch ${VERSION} --single-branch https://github.com/luntergroup/octopus.git \
    && cd octopus \
    # Patch to avoid compiling 
    # with AVX2 instruction set
    # when building the image
    # on a modern machine, this
    # instruction set will fail 
    # on nodes that have older 
    # generation CPUs, ivybridge
    # or sandybridge, here is a 
    # list of CPUs that support 
    # this instuction set:
    # https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#Advanced_Vector_Extensions_2 
    && sed -i 's/^elseif (AVX2_FOUND)/#elseif (AVX2_FOUND)/g' src/CMakeLists.txt \
    && sed -i 's/add_compile_options(-mavx2)/#add_compile_options(-mavx2)/g' src/CMakeLists.txt \
    && mkdir -pv build_dir \
    && cd build_dir \
    && cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCOMPILER_ARCHITECTURE=sandybridge \
        -DHTSLIB_ROOT=/usr/lib/x86_64-linux-gnu \
        .. \
    && make install -j 1

# Add Dockerfile and export env variables
WORKDIR /data2
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2/octopus/bin:${PATH}"
WORKDIR /data2