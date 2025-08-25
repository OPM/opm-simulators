#!/bin/bash
set -eux

# These should be passed as environment variables from the Dockerfile.
: "${OPM_COMMON_REPO:?}"
: "${OPM_SIMULATORS_REPO:?}"
: "${OPM_COMMON_BRANCH:?}"
: "${OPM_SIMULATORS_BRANCH:?}"

git clone --depth=1 "$OPM_COMMON_REPO" -b "$OPM_COMMON_BRANCH" opm-common
git clone --depth=1 "$OPM_SIMULATORS_REPO" -b "$OPM_SIMULATORS_BRANCH" opm-simulators

mv opm-simulators/python/opm opm-simulators/python/opm_original
mv opm-common/python/opm opm-common/python/opm_original
