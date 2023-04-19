############# OVERVIEW #############
# genome-seek sv docker image,
# uses genome-seek as a base image:
FROM skchronicles/genome-seek:v0.1.0

# Base image already contains:
# bzip2, gcc, g++, make, python, zlib1g-dev 
# @INSTALLS:
#  • manta/1.6.0
WORKDIR /opt2


############## INSTALL #############
# Install manta/1.6.0, for more info see:
# https://github.com/Illumina/manta
RUN mkdir -p /opt2/manta/1.6.0/ \
    && wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 -P /opt2/manta/1.6.0/ \
    && tar -xvjf /opt2/manta/1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 -C /opt2/manta/1.6.0/ \
    && rm /opt2/manta/1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && rm -rf /opt2/manta/1.6.0/manta-1.6.0.centos6_x86_64/share/demo/ \
    && ln -fs /usr/bin/python2.7 /usr/bin/python2
ENV PATH="${PATH}:/opt2/manta/1.6.0/manta-1.6.0.centos6_x86_64/bin"


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek_sv.dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2
