#!/bin/bash
set -eux

# This script is called from the Dockerfiles, e.g. docker/test_wheels/dockerfiles/.../Dockerfile
# It is used to install and setup pyenv

# NOTE: This variable is set in the Dockerfile so we do not need to set it here.
# export PYENV_ROOT=/root/.pyenv
: "${PYENV_ROOT:?}"

git clone https://github.com/pyenv/pyenv.git "$PYENV_ROOT"

echo 'export PYENV_ROOT="/root/.pyenv"' >> /root/.bashrc
echo 'export PATH="${PYENV_ROOT}/bin:${PYENV_ROOT}/shims:${PATH}"' >> /root/.bashrc
echo 'eval "$(pyenv init - bash)"' >> /root/.bashrc
