#!/usr/bin/env bash
set -e

pushd . > /dev/null
cd opm-material
mkdir build
cd build
cmake ../
make
popd > /dev/null
