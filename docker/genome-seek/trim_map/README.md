## Steps for Building Docker Images

This docker image uses the following [base image](https://github.com/OpenOmics/genome-seek/blob/main/docker/genome-seek/Dockerfile).

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=genome-seek_trim_map:v0.1.1 .

# Testing, take a peek inside
docker run -ti genome-seek_trim_map:v0.1.1 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag genome-seek_trim_map:v0.1.1 skchronicles/genome-seek_trim_map:v0.1.1
docker tag genome-seek_trim_map:v0.1.1 skchronicles/genome-seek_trim_map         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/genome-seek_trim_map:v0.1.1
docker push skchronicles/genome-seek_trim_map:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan genome-seek_trim_map:v0.1.1
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
