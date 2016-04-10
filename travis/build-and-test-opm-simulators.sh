#!/usr/bin/env bash
set -ex

pushd . > /dev/null
opm-simulators/travis/build-opm-simulators.sh
cd opm-simulators/build
ctest --output-on-failure
popd > /dev/null
