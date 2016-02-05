#!/usr/bin/env bash
set -e

pushd . > /dev/null
git clone https://github.com/dune-project/dune-common.git
cd dune-common
git checkout tags/v2.3.1
mkdir build
cd build
cmake ../
make
popd > /dev/null
