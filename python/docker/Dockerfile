# Dockerfile to generate PyPI packages. Needs to be run from the opm-simulators/python folder
#
# NOTE: Preferred way to build this docker image is to use docker/run-docker-build.sh
#
# Example use:
# docker build -t manylinux2014_opm:static . -f docker/Dockerfile
# docker build -t manylinux2014_opm:shared . -f docker/Dockerfile --build-arg libtype=shared
# docker build -t manylinux2014_opm:py312 . -f docker/Dockerfile --build-arg python_versions="3.12"
# docker build -t manylinux2014_opm:pr-test . -f docker/Dockerfile --build-arg version_simulators="pull/1234/head"
#
# To extract the wheel files from the final scratch image to a local directory, add --output:
# docker build -t manylinux2014_opm:static . -f docker/Dockerfile --output wheelhouse
# docker build -t manylinux2014_opm:shared . -f docker/Dockerfile --build-arg libtype=shared --output wheelhouse-shared
#
# Note: To see detailed build output instead of buildx progress, use:
# DOCKER_BUILDKIT=0 docker build -t manylinux2014_opm:built . -f docker/Dockerfile
# To redirect output to a log file for later inspection:
# docker build -t manylinux2014_opm:shared . -f docker/Dockerfile --build-arg libtype=shared 2>&1 | tee build.log

# Using manylinux_2_28 for broad compatibility (Debian 10+, Ubuntu 18.10+, CentOS/RHEL 8+)
# Newer manylinux_2_34 is available but has narrower distribution support
FROM quay.io/pypa/manylinux_2_28_x86_64 AS stage1
# Version control for each OPM repository (supports branches, tags, or PR refs like "pull/1234/head")
ARG version_common="master"
ARG version_grid="master"
ARG version_simulators="master"
ARG libtype="static"  # "static" for static libraries (.a), "shared" for shared libraries (.so)
WORKDIR /tmp/opm
RUN if [ "$libtype" != "static" ] && [ "$libtype" != "shared" ]; then \
        echo "Error: libtype must be 'static' or 'shared', got '$libtype'"; \
        exit 1; \
    fi
RUN echo "Using versions - Common: $version_common, Grid: $version_grid, Simulators: $version_simulators"
# Copy only the setup script to avoid invalidating cache with source code changes
COPY docker/scripts/setup-docker-image.sh /tmp/opm/setup-docker-image.sh
RUN /bin/bash setup-docker-image.sh $libtype

# Multi-stage build: stage1=setup, stage2=build deps, stage3=generate packages, export_stage=extract results
FROM stage1 AS stage2
ARG libtype="static"
ARG build_jobs=4
# Build dependencies in separate layers for better Docker caching
# Copy only dependency build scripts (delay "ADD . .") to preserve cache layers for better debugging
#   experience. I.e. if you change the source code, the build will not be invalidated and you have to
#   rebuild everything from scratch.
COPY docker/scripts/build-boost.sh /tmp/opm/build-boost.sh
RUN /bin/bash build-boost.sh $build_jobs $libtype
COPY docker/scripts/build-dune.sh /tmp/opm/build-dune.sh
RUN /bin/bash build-dune.sh $build_jobs $libtype

FROM stage2 AS stage3
# Version control for each OPM repository
ARG version_common="master"
ARG version_grid="master"
ARG version_simulators="master"
ARG libtype="static"
ARG build_jobs=4
ARG version_tag=""
# CMake build targets (default for master branch, may differ in PRs)
ARG target_common="opmcommon_python"
ARG target_simulators="simulators"
# Default Python versions - can be overridden via build-arg or read from docker/python_versions.json
# NOTE: The list is maintained in docker/python_versions.json and docker/scripts/sync_versions.sh
ARG python_versions="3.8,3.9,3.10,3.11,3.12,3.13"
# Finally copy all source files - only needed for wheel generation
ADD . .
RUN /bin/bash docker/scripts/generate-pypi-package.sh "$version_common" "$version_grid" "$version_simulators" "$build_jobs" "$version_tag" $libtype "$python_versions" "$target_common" "$target_simulators"

# Final stage: extract only the wheel files using minimal scratch image
FROM scratch AS export_stage
COPY --from=stage3 /tmp/opm/wheelhouse .
