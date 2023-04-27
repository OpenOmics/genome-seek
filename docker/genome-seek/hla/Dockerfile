############# OVERVIEW #############
# genome-seek hla docker image,
# uses genome-seek as a base image:
FROM skchronicles/genome-seek:v0.1.0
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="skchronicles/genome-seek:v0.1.0" \
    version="v0.1.0"   \
    software="genome-seekq hla" \
    about.summary="An awesome Clinical Whole Genome Sequencing Pipeline" \
    about.home="https://github.com/OpenOmics/genome-seek" \
    about.documentation="https://openomics.github.io/genome-seek/" \
    about.tags="Human leukocyte antigen (HLA) type inference"

# Base image already contains:
# g++ (support for C++11), libboost (apt-get version not compatible with python3.10), 
# zlib, samtools (1.17), picard
# @INSTALLS:
#  • bamtools/2.5.2  ✓
#  • bwa/0.7.17-6    ✓
#  • libboost/1.80.0 ✓
#  • HLA-LA/1.0.3    ✓
#  • bamtools/2.5.2  ✓
WORKDIR /opt2
# Install bwa/0.7.17-6 for HLA-LA,
# perl modules for HLA-LA,
# libicu-dev for boost
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        autotools-dev \
        bwa \
        bioperl \
        libicu-dev \
        libbio-db-hts-perl \
        libbio-featureio-perl \
        libbio-samtools-perl \
        libdata-dump-perl \
        libjsoncpp25 \
        liblist-moreutils-perl \
        libmath-basecalc-perl \
        libparallel-forkmanager-perl \
        libtext-levenshtein-perl \
    && cpan Array::Utils \
    && perl -MArray::Utils -e ';' \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


############## INSTALL #############
# Install C++ Boost/1.80.0 library,
# required for HLA-LA, apt install
# of libboost libraries are missing
# include/ and install/ directories 
# needed to make the HLA-LA exe, 
# also, versions of boost<1.80.0 
# are NOT compatiable with python3.10,
# for more information see: 
# https://www.boost.org/doc/libs/1_80_0/more/getting_started/unix-variants.html
# http://www.boost.org/build/doc/html/index.html
# https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=1000506
# https://github.com/boostorg/python/pull/344
RUN mkdir -p /opt2/lib/ \
    && wget "http://downloads.sourceforge.net/project/boost/boost/1.80.0/boost_1_80_0.tar.bz2" -P /opt2/lib/ \
    && tar -xjf /opt2/lib/boost_1_80_0.tar.bz2 -C /opt2/lib/ \
    && rm /opt2/lib/boost_1_80_0.tar.bz2 \
    && cd /opt2/lib/boost_1_80_0 \
    && ./bootstrap.sh \
    && ./b2 -q install --prefix=/opt2/libboost 
ENV PATH="${PATH}:/opt2/libboost/bin" \
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt2/libboost/lib"

# Install bamtools/2.5.2,
# for more info see:
# https://github.com/pezmaster31/bamtools
# https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
RUN mkdir -p /opt2/bamtools/ \
    && wget https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.5.2.tar.gz -P /opt2/bamtools/ \
    && tar -vxzf /opt2/bamtools/v2.5.2.tar.gz -C /opt2/bamtools/ \
    && rm /opt2/bamtools/v2.5.2.tar.gz \
    && mkdir -p /opt2/bamtools/bamtools-2.5.2/build \
    && cd /opt2/bamtools/bamtools-2.5.2/build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt2/libbamtools \
        -DCMAKE_INSTALL_INCLUDEDIR=/opt2/libbamtools/include \
        -DCMAKE_INSTALL_LIBDIR=/opt2/libbamtools/lib \
        -DBUILD_SHARED_LIBS=TRUE \
        .. \
    && make install \
    && rm -rf /opt2/bamtools/bamtools-2.5.2/build
ENV PATH="${PATH}:/opt2/libbamtools/bin" \
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt2/libbamtools/lib"

# Install HLA-LA/1.0.3,
# for more info see:
# https://github.com/DiltheyLab/HLA-LA
# https://github.com/DiltheyLab/HLA-LA#installing-hlala
# Checkout commit tied to v1.0.3
# and test that the binary works!
RUN mkdir -p /opt2/hla-la/1.0.3/ \
        /opt2/hla-la/1.0.3/HLA-LA \
        /opt2/hla-la/1.0.3/HLA-LA/bin \
        /opt2/hla-la/1.0.3/HLA-LA/src \
        /opt2/hla-la/1.0.3/HLA-LA/obj \
        /opt2/hla-la/1.0.3/HLA-LA/temp \
        /opt2/hla-la/1.0.3/HLA-LA/working \
        /opt2/hla-la/1.0.3/HLA-LA/graphs \
    && cd /opt2/hla-la/1.0.3/HLA-LA/src \
    && git clone https://github.com/DiltheyLab/HLA-LA.git . \
    && git checkout fb8ae1cde7783d6cf441a80319e943f63f3b6276 \
    && ln -sf /opt2/libbamtools/lib /opt2/libbamtools/lib64 \
    && make all \
        BOOST_PATH=/opt2/libboost \
        BAMTOOLS_PATH=/opt2/libbamtools \
    && sed -i "s:^\(picard_sam2fastq_bin=\).*:\1$PICARD:" paths.ini \
    && sed -i "s:^\(workingDir=\).*:workingDir=:" paths.ini \
    && sed -i "s:^\(workingDir_HLA_ASM=\).*:workingDir_HLA_ASM=:" paths.ini \
    && /opt2/hla-la/1.0.3/HLA-LA/bin/HLA-LA --action testBinary
ENV PATH="${PATH}:/opt2/hla-la/1.0.3/HLA-LA/bin:/opt2/hla-la/1.0.3/HLA-LA/src"


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek_hla.dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2