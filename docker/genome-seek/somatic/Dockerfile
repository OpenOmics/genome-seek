############# OVERVIEW #############
# genome-seek somatic docker image,
# uses genome-seek as a base image:
FROM skchronicles/genome-seek:v0.1.0
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="skchronicles/genome-seek:v0.1.0" \
    version="v0.1.0"   \
    software="genome-seekq somatic" \
    about.summary="An awesome Clinical Whole Genome Sequencing Pipeline" \
    about.home="https://github.com/OpenOmics/genome-seek" \
    about.documentation="https://openomics.github.io/genome-seek/" \
    about.tags="Variant-calling quality-control"

# Base image already contains:
# gatk3 (add GATK wrapper from HPC staff),
# gatk4, bcftools, vcftools, gatk4,
# R/4.X, R-maftools
# Octopus is in its own seperate image.
# @INSTALLS:
#  • muse       ✓
#  • strelka    ✓ (indices include in resource bundle)
WORKDIR /opt2


############## INSTALL #############
# Add GATK wrapper to gatk3, provided
# by Biowulf HPC staff. gatk3 is 
# installed in the base image.
RUN mkdir -p /opt2/gatk3/
COPY src/GATK /opt2/gatk3/
ENV PATH="${PATH}:/opt2/gatk3"

# Install muse/2.0.1, for more info see:
# https://github.com/wwylab/MuSE/tree/master
# Checking out exact commit that is associated
# with the v2.0.1 release, not sure if release
# tarball contains everything.
RUN mkdir -p /opt2/muse/2.0.1/ \
    && git clone --recursive https://github.com/wwylab/MuSE /opt2/muse/2.0.1/MuSE \
    && cd /opt2/muse/2.0.1/MuSE \
    && git checkout 105c5be01c43dbe8cc29803e2ab87f35f6f3589a \
    && ./install_muse.sh \
    && mkdir bin \
    && mv /opt2/muse/2.0.1/MuSE/MuSE /opt2/muse/2.0.1/MuSE/bin
ENV PATH="${PATH}:/opt2/muse/2.0.1/MuSE/bin"
WORKDIR /opt2

# Install strelka/2.9.10, for more info see:
# https://github.com/Illumina/strelka
# Set python2.7 as the default interpreter
# @Required system libraries:
#  • python 2.6+     ✓ Already in base image
#  • bzip2           ✓ Already in base image
#  • gcc             ✓ Already in base image 
#  • g++             ✓ Already in base image 
#  • make            ✓ Already in base image 
#  • zlib1g-dev      ✓ Already in base image
RUN mkdir -p /opt2/strelka/2.9.10/ \
    && wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 -P /opt2/strelka/2.9.10/ \
    && tar -xvjf /opt2/strelka/2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 -C /opt2/strelka/2.9.10/ \
    && rm /opt2/strelka/2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && ln -fs /usr/bin/python2.7 /usr/bin/python2
ENV PATH="${PATH}:/opt2/strelka/2.9.10/strelka-2.9.10.centos6_x86_64/bin"


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek_somatic.dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2