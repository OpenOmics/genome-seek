############# OVERVIEW #############
# genome-seek trim_map docker image,
# uses previous version of genome-seek 
# trim_map (v0.1.0) as a base image:
FROM skchronicles/genome-seek_trim_map:v0.1.0

# Base image already contains:
# samblaster, samtools, fastp, bwa_mem2
WORKDIR /opt2

############## INSTALL ##############
# Add bwa-mem2-wrapper to select the
# correct bwa-mem2 pre-compiled binary
# for Intel and AMD (Zen4) processors.
# AMD (Zen4) processors require are not
# compatible with the AVX-512bw Intel 
# pre-compiled binary. This wrapper will
# select the AVX2 binary if run on an
# AMD processor; otherwise, it will
# let bwa-mem2 select the best binary.
ADD bin/bwa-mem2-wrapper /opt2/bwa_mem2/2.2.1/bwa-mem2-2.2.1_x64-linux/bwa-mem2-wrapper
RUN chmod a+rx /opt2/bwa_mem2/2.2.1/bwa-mem2-2.2.1_x64-linux/bwa-mem2-wrapper

################ POST #################
# Add Dockerfile and export env variables
ADD Dockerfile /opt2/genome-seek_trim_map_v0.1.1.dockerfile
WORKDIR /data2