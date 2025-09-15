#!/bin/bash

# This script is called from docker/Dockerfile to build the boost library.

set -e

BUILD_JOBS=$1
LIBTYPE=$2

if [ "$LIBTYPE" == "static" ]; then
    shared=0
else
    shared=1
fi

export CMAKE_GENERATOR=Ninja

pushd /tmp/opm

# Build boost
git clone --depth 1 --branch boost-1.84.0 https://github.com/boostorg/boost
pushd boost
git submodule init
git submodule update
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=$shared -DCMAKE_POSITION_INDEPENDENT_CODE=1 -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=1
cmake --build . -j ${BUILD_JOBS}
cmake --build . --target install
popd
popd
