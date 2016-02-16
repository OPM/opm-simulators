#!/usr/bin/env bash
set -ex

pushd . > /dev/null
cd opm-autodiff
mkdir build
cd build
cmake -D SUPERLU_ROOT=../../SuperLU ../
make
popd > /dev/null
