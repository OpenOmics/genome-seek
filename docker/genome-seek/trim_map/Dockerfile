############# OVERVIEW #############
# genome-seek trim_map docker image,
# uses genome-seek as a base image:
FROM skchronicles/genome-seek:v0.1.0

# Base image already contains:
# samblaster and samtools
# @INSTALLS: 
#  • fastp, bwa_mem2
WORKDIR /opt2

############## INSTALL #############
# Install fastp, for more info see:
# https://github.com/OpenGene/fastp/
RUN mkdir -p /opt2/fastp/0.23.2/ \
    && wget http://opengene.org/fastp/fastp.0.23.2 -P /opt2/fastp/0.23.2/ \
    && mv /opt2/fastp/0.23.2/fastp.0.23.2 /opt2/fastp/0.23.2/fastp \
    && chmod -R a+rx /opt2/fastp/0.23.2/fastp
ENV PATH="${PATH}:/opt2/fastp/0.23.2"

# Install bwa_mem2, for more info see:
# https://github.com/bwa-mem2/bwa-mem2
RUN mkdir -p /opt2/bwa_mem2/2.2.1/ \
    && wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 -P /opt2/bwa_mem2/2.2.1/ \
    && tar -xjvf /opt2/bwa_mem2/2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 -C /opt2/bwa_mem2/2.2.1/ \
    && rm /opt2/bwa_mem2/2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
ENV PATH="${PATH}:/opt2/bwa_mem2/2.2.1/bwa-mem2-2.2.1_x64-linux"


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek_trim_map.dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2
