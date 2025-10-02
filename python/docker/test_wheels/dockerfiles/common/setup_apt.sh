#! /bin/bash

set -eux

# This script is called from the Dockerfiles, e.g. docker/test_wheels/dockerfiles/.../Dockerfile
# It is used to setup the apt repositories and install the necessary packages for debian based images.
#
# NOTE: These variables are set in the Dockerfile so we do not need to set them here.
#export DEBIAN_FRONTEND=noninteractive
#export TZ=${TZ:-Etc/UTC}
: "${DEBIAN_FRONTEND:?}"
: "${TZ:?}"

apt-get update
apt-get install -y --no-install-recommends apt-utils tzdata

ln -snf /usr/share/zoneinfo/$TZ /etc/localtime
echo $TZ > /etc/timezone

apt-get install -y --no-install-recommends \
    python3 python3-pip python3-dev \
    git build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev \
    libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils \
    tk-dev libffi-dev liblzma-dev

rm -rf /var/lib/apt/lists/*
