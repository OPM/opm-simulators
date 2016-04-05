#!/bin/bash

# Build flow binary
pushd .
cd serial/build-opm-autodiff
cmake --build . --target flow
popd

# Clone opm-data if necessary
pushd .
cd deps
if ! test -d opm-data
then
  git clone --depth 1 --single-branch -b master https://github.com/OPM/opm-data
fi
cd opm-data

# Run the norne case
cd norne
$WORKSPACE/serial/build-opm-autodiff/bin/flow deck_filename=NORNE_ATW2013.DATA output_dir=OPM
test $? -eq 0 || exit 1
PATH=$WORKSPACE/serial/install/bin:$PATH ./plotwells.sh

popd
