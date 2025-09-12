FROM debian:11

# Set build-time and runtime environment
# Avoid user-interaction prompts when installing packages.
ARG DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_FRONTEND=${DEBIAN_FRONTEND}
ENV TZ=Etc/UTC
ENV PYENV_ROOT=/root/.pyenv
ENV PATH="${PYENV_ROOT}/bin:${PYENV_ROOT}/shims:${PATH}"

WORKDIR /test

# Make all scripts in the common folder available
COPY test_wheels/dockerfiles/common/ /test/common/
# Copy python_versions.json to the expected location
COPY python_versions.json /test/python_versions.json

RUN /test/common/setup_apt.sh
RUN /test/common/setup_pyenv.sh

# Switch to bash for pyenv usage
SHELL ["/bin/bash", "-c"]

RUN /usr/bin/python3 /test/common/install_python_versions.py

# Clone opm-common and opm-simulators to get the unit tests
ARG OPM_COMMON_REPO=https://github.com/OPM/opm-common.git
ARG OPM_SIMULATORS_REPO=https://github.com/OPM/opm-simulators.git
ARG OPM_COMMON_BRANCH=master
ARG OPM_SIMULATORS_BRANCH=master
ENV OPM_COMMON_REPO=${OPM_COMMON_REPO}
ENV OPM_SIMULATORS_REPO=${OPM_SIMULATORS_REPO}
ENV OPM_COMMON_BRANCH=${OPM_COMMON_BRANCH}
ENV OPM_SIMULATORS_BRANCH=${OPM_SIMULATORS_BRANCH}
RUN /test/common/clone_opm_repos.sh

# Entrypoint: use unbuffered mode so logs appear in order
ENTRYPOINT ["/usr/bin/python3", "-u", "/test/common/entrypoint.py"]
# No default arguments passed to entrypoint
CMD []
