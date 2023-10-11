## Steps for Building Docker Images

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=superfreq:v0.1.1 .

# Testing, take a peek inside
docker run -ti superfreq:v0.1.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag superfreq:v0.1.1 skchronicles/superfreq:v0.1.1
docker tag superfreq:v0.1.1 skchronicles/superfreq         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/superfreq:v0.1.1
docker push skchronicles/superfreq:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan superfreq:v0.1.1
```

### Gotchas

[VEP's reference cache and plugins](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html) are not installed with the container's filesystem. We recommend using a local cache to avoid making any external requests(API may not be reliable). As so, please download your reference cache of interest and any plugins you may require prior to run time (don't forget to bind them to the container's filesystem). superFreq provides a wrapper to VEP to annotate the results. The [runVEP.R module provided superFreq](https://github.com/ChristofferFlensburg/superFreq/blob/master/R/runVEP.R) does not have any built-in options to specify the location of the cache/plugins directory, so it defaults to `~/.vep/`. This is problematic unless your cache directory is in this location.

To avoid this problem, please export the following environment variables prior to running superFreq:
```bash
# Important/required VEP enviroment variables
export VEP_DIR_CACHE=/path/to/VEP/110/cache/
export VEP_OFFLINE=1
export VEP_FORCE_OVERWRITE=1
export VEP_CACHE=1

# Setting up tmpdir for Rscripts,
# tempdir() function in R will use
# these variables to resolve the 
# location writing temp files
mkdir -p /path/to/superFreq/outdir
export TMPDIR=/path/to/superFreq/outdir
export TMP=/path/to/superFreq/outdir
export TEMP=/path/to/superFreq/outdir
```

This will ensure [superFreq](https://github.com/ChristofferFlensburg/superFreq) runs `VEP/110` in an offline mode using your local cache. 

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
