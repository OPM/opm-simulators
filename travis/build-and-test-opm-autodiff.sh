#!/usr/bin/env bash
set -ex

pushd . > /dev/null
opm-autodiff/travis/build-opm-autodiff.sh
cd opm-autodiff/build
ctest --output-on-failure
popd > /dev/null
