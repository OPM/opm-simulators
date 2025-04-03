#!/bin/bash

# Script to be run on a manylinux2014 docker image to complete it for OPM usage.
# i.e. docker run -i -t quay.io/pypa/manylinux2014_x86_64 < setup-docker-image.sh

# A ready made Docker image is available at Dockerhub:
# docker run -i -t lindkvis/manylinux2014_opm:latest

LIBTYPE=${1:-"static"}

dnf install -y almalinux-release-devel ninja-build

if [ "$LIBTYPE" == "static" ]; then
    dnf install -y blas-static lapack-static suitesparse-static
else
    dnf install -y blas-devel lapack-devel suitesparse-devel
fi
