# Base image using python/3.8
# Avoids issue described here:
# https://github.com/KarchinLab/open-cravat/issues/113
FROM python:3.8

LABEL maintainer=kuhnsa@nih.gov
ARG PACKAGEDIR=/usr/local/lib/python3.6/site-packages/cravat
ARG VERSION=2.2.7

# Create Container filesystem specific 
# working directory and opt directories 
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
        git \
        locales \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Install OpenCRAVAT from official release 
RUN wget https://github.com/KarchinLab/open-cravat/archive/refs/tags/${VERSION}.tar.gz \
    && tar -xvf ${VERSION}.tar.gz \
    && rm ${VERSION}.tar.gz \
    && cd open-cravat-${VERSION} \
    && pip install . \
    && pip install open-cravat-multiuser aiosqlite3 scipy \
    && oc version

WORKDIR /opt2 

# Add Dockerfile and export env variables
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2 
ENV PATH="/opt2:$PATH"
WORKDIR /data2
