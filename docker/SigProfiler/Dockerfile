# SigProfiler Dockerfile
# Using Ubuntu Jammy (22.04 LTS) as base image.
# SigProfiler is a suite of seperate tools with 
# that perform different functions:
# https://cancer.sanger.ac.uk/signatures/tools/
# This docker image contains:
#  • SigProfilerMatrixGenerator:
#    • https://github.com/AlexandrovLab/SigProfilerMatrixGenerator
#  • SigProfilerExtractor:
#    • https://github.com/AlexandrovLab/SigProfilerExtractor
#  • SigProfilerSimulator:
#    • https://github.com/AlexandrovLab/SigProfilerSimulator
#  • SigProfilerPlotting:
#    • https://github.com/AlexandrovLab/SigProfilerPlotting
FROM ubuntu:22.04
LABEL maintainer=kuhnsa@nih.gov

# Overview of Dependencies
#  • python version >= 3:   python3 (3.10.4)
#  • Python3 Package(s):
#    • SigProfilerMatrixGenerator
#    • SigProfilerExtractor
#    • SigProfilerSimulator
#    • SigProfilerPlotting

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
        curl \
        git \
        locales \
        python3 \
        python3-pip \
        rsync \
        wget \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Install Python Packages
RUN ln -s /usr/bin/python3.10 /usr/bin/python
RUN pip3 install --upgrade pip \
	  && pip3 install numpy pandas matplotlib argparse \
      && pip3 install \
            SigProfilerMatrixGenerator \
            SigProfilerExtractor \
            SigProfilerSimulator \
            SigProfilerPlotting

# Install GRCh38 reference files,
# In future more refs maybe addded.
RUN wget http://ngs.sanger.ac.uk/scratch/project/mutographs/SigProf/GRCh38.tar.gz \
    && echo 'from SigProfilerMatrixGenerator import install as genInstall' > /opt2/install.py \
    && echo 'genInstall.install("GRCh38", offline_files_path="/opt2/")' >> /opt2/install.py \
    && python3 /opt2/install.py \
    && rm -f /opt2/GRCh38.tar.gz

# Add Dockerfile and export env variables
WORKDIR /data2
ADD Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2