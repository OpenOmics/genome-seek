############# OVERVIEW #############
# genome-seek qc docker image,
# uses genome-seek as a base image:
FROM skchronicles/genome-seek:v0.1.0
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>" \
    base_image="skchronicles/genome-seek:v0.1.0" \
    version="v0.1.0"   \
    software="genome-seekq qc" \
    about.summary="An awesome Clinical Whole Genome Sequencing Pipeline" \
    about.home="https://github.com/OpenOmics/genome-seek" \
    about.documentation="https://openomics.github.io/genome-seek/" \
    about.tags="Variant-calling Quality-control"

# Base image already contains:
# bcftools, vcftools, python, perl, samtools, 
# java 8, gatk4, picard, bowtie2/2.4.4, R/4.X
# @INSTALLS:
#  • fastqc       ✓
#  • fastq_screen ✓ (indices include in resource bundle)
#  • qualimap     ✓
#  • snpeff       ✓
#  • somalier     ✓
#  • multiqc/1.14 ✓
WORKDIR /opt2
# Adding depedencies for fastq_screen GD::Graph::bars
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        libgd-perl \
        pkg-config \
        libgd-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

############## INSTALL #############
# Install MultiQC/1.14, for more info see:
RUN pip3 install --upgrade pip \
	  && pip3 install multiqc==1.14

# Install fastqc, for more info see:
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
RUN mkdir -p /opt2/fastqc/0.12.1/ \
    && cd /opt2/fastqc/0.12.1/ \
    && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip -P /opt2/fastqc/0.12.1/ \
    && unzip /opt2/fastqc/0.12.1/fastqc_v0.12.1.zip \
    && rm /opt2/fastqc/0.12.1/fastqc_v0.12.1.zip
ENV PATH="${PATH}:/opt2/fastqc/0.12.1/FastQC"
WORKDIR /opt2

# Install fastq_screen, for more info see:
# https://stevenwingett.github.io/FastQ-Screen/
RUN mkdir -p /opt2/fastq_screen/0.15.3/ \
    && cpanm GD::Graph::bars \
    && wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz -P /opt2/fastq_screen/0.15.3/ \
    && tar -vxzf /opt2/fastq_screen/0.15.3/v0.15.3.tar.gz -C /opt2/fastq_screen/0.15.3/ \
    && rm /opt2/fastq_screen/0.15.3/v0.15.3.tar.gz
ENV PATH="${PATH}:/opt2/fastq_screen/0.15.3/FastQ-Screen-0.15.3"

# Install Qualimap/2.2.1, for more info see:
# http://qualimap.conesalab.org/
# https://gensoft.pasteur.fr/docs/qualimap/v2.2.1/QualimapManual.pdf
# @Required R packages:
#  • optparse (CRAN)                ✓ Already in base image
#  • NOISeq  (Bioconductor)         - Needs to be installed
#  • Repitools (Bioconductor)       - Needs to be installed
#  • Rsamtools (Bioconductor)       ✓ already in base image
#  • GenomicFeatures (Bioconductor) ✓ already in base image
#  • rtracklayer (Bioconductor)     ✓ already in base image
RUN mkdir -p /opt2/qualimap/2.2.1 \
    && Rscript -e "BiocManager::install(c('NOISeq', 'Repitools'), Ncpus = 4)" \
    && wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip -P /opt2/qualimap/2.2.1/ \
    && cd /opt2/qualimap/2.2.1/ \
    && unzip /opt2/qualimap/2.2.1/qualimap_v2.2.1.zip \
    && rm /opt2/qualimap/2.2.1/qualimap_v2.2.1.zip
ENV PATH="${PATH}:/opt2/qualimap/2.2.1/qualimap_v2.2.1"
WORKDIR /data2

# Install somalier/0.2.16, for more info see:
# https://github.com/brentp/somalier
RUN mkdir -p /opt2/somalier/0.2.16 \
    && wget https://github.com/brentp/somalier/releases/download/v0.2.16/somalier -P /opt2/somalier/0.2.16/ \
    && chmod a+rx /opt2/somalier/0.2.16/somalier
ENV PATH="${PATH}:/opt2/somalier/0.2.16"

# Install snpEff/4.3t, for more info see:
# https://pcingola.github.io/SnpEff/
# Version 4.3t works with java 8 (already in base image),
# deleting bundled version of clinEff (not being used),
# setting env variable to jar file similar to HPC module
RUN mkdir -p /opt2/snpEff/4.3t \
    && cd /opt2/snpEff/4.3t/ \
    && wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip -P /opt2/snpEff/4.3t/ \
    && unzip /opt2/snpEff/4.3t/snpEff_v4_3t_core.zip \
    && rm /opt2/snpEff/4.3t/snpEff_v4_3t_core.zip \
    && rm -rf /opt2/snpEff/4.3t/clinEff/
ENV SNPEFF_JAR="/opt2/snpEff/4.3t/snpEff/snpEff.jar"


################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek_qc.dockerfile
RUN chmod -R a+rX /opt2
WORKDIR /data2