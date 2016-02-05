#!/usr/bin/env bash
set -e

pushd . > /dev/null
opm-material/travis/build-opm-material.sh
cd opm-material/build
ctest --output-on-failure
popd > /dev/null
