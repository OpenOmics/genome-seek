## Steps for Building Docker Images

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f example.dockerfile --tag=example:v0.1.0 .

# Testing, take a peek inside
docker run -ti example:v0.1.0 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag example:v0.1.0 skchronicles/example:v0.1.0
docker tag example:v0.1.0 skchronicles/example         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/example:v0.1.0
docker push skchronicles/example:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan example:v0.1.0
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
