FROM ubuntu:22.04

# Avoid user-interaction prompts when installing packages.
ARG DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_FRONTEND=${DEBIAN_FRONTEND}
ENV TZ=Etc/UTC

# Install system dependencies for building Python from source
RUN apt-get update && apt-get install -y \
 && apt-get install -y apt-utils \
 && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
 && echo $TZ > /etc/timezone \
 && apt-get install -y \
    python3 python3-pip python3-dev \
    git build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev \
    libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils \
    tk-dev libffi-dev liblzma-dev \
 && rm -rf /var/lib/apt/lists/*

# Set up pyenv
ENV PYENV_ROOT /root/.pyenv
# Put both 'bin' and 'shims' on PATH for everything
ENV PATH="${PYENV_ROOT}/bin:${PYENV_ROOT}/shims:${PATH}"

RUN git clone https://github.com/pyenv/pyenv.git $PYENV_ROOT

# We’ll explicitly use bash for the following commands
SHELL ["/bin/bash", "-c"]

# Load pyenv automatically in future shells
RUN echo 'export PYENV_ROOT="${PYENV_ROOT}"' >> /root/.bashrc \
 && echo 'export PATH="${PYENV_ROOT}/bin:${PYENV_ROOT}/shims:${PATH}"' >> /root/.bashrc \
 && echo 'eval "$(pyenv init - bash)"' >> /root/.bashrc

# Map from major python version to the specific patch version to be installed by pyenv
ENV PYTHON_VERSION_MAP='{"3.6": "3.6.15", "3.7": "3.7.17", "3.8": "3.8.19", "3.9": "3.9.21", "3.10": "3.10.16", "3.11": "3.11.11", "3.12": "3.12.9"}'

COPY dockerfiles/ubuntu:22.04/install_python_versions.py .
RUN /usr/bin/python3 install_python_versions.py

WORKDIR /test

# Clone opm-common and opm-simulators to get the unit tests
ARG opm_common_repo=https://github.com/OPM/opm-common.git
ARG opm_simulators_repo=https://github.com/OPM/opm-simulators.git
ARG opm_common_branch=master
ARG opm_simulators_branch=master
RUN git clone --depth=1 $opm_common_repo -b $opm_common_branch opm-common
RUN git clone --depth=1 $opm_simulators_repo -b $opm_simulators_branch opm-simulators
# Rename the python/opm folders to avoid confusion the globally installed opm package
#  from the wheels
RUN mv opm-simulators/python/opm opm-simulators/python/opm_original
RUN mv opm-common/python/opm opm-common/python/opm_original

# The wheel directory will be mounted at runtime at /test/wheels

# Copy in the custom entrypoint script
COPY dockerfiles/ubuntu:22.04/entrypoint.py .
RUN chmod +x entrypoint.py

# Run our entrypoint on container start
# Use "-u" to force stdout/stderr to be unbuffered. This is important for the python
# script output to be visible in terminal immediately, else the output from python's
# subprocess.run() will always appear before any print() output and the joined output will be
# unreadable.
ENTRYPOINT ["/usr/bin/python3", "-u", "/test/entrypoint.py"]
CMD []
