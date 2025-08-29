#!/bin/bash

# This script is used to setup a manylinux docker image for OPM usage. It is called from
#   docker/Dockerfile

LIBTYPE=${1:-"static"}

dnf install -y almalinux-release-devel ninja-build

if [ "$LIBTYPE" == "static" ]; then
    dnf install -y blas-static lapack-static suitesparse-static
else
    dnf install -y blas-devel lapack-devel suitesparse-devel
fi
