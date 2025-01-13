# Dockerfile to generate PyPI packages. Needs to be run from the opm-common root folder
# Example use:
# sudo docker build -t manylinux2014_opm:built . -f python/Dockerfile

FROM quay.io/pypa/manylinux_2_28_x86_64 AS stage1
ARG version="master"
ARG build_jobs=4
ARG version_tag=""
WORKDIR /tmp/opm
RUN echo "Using version: $version"
RUN echo "Using package version tag: $version_tag"
RUN echo "Number of build jobs: $build_jobs"
ADD . .
RUN /bin/bash setup-docker-image.sh
FROM stage1 AS stage2
RUN /bin/bash build-deps.sh $build_jobs
FROM stage2 AS stage3
RUN /bin/bash generate-pypi-package.sh "$version" "$build_jobs" "$version_tag"

FROM scratch AS export_stage
COPY --from=stage3 /tmp/opm/wheelhouse .
