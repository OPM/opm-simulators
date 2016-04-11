#!/bin/bash

pushd .
cd deps/opm-data

# Run the norne case
cd norne
$WORKSPACE/serial/build-opm-simulators/bin/flow deck_filename=NORNE_ATW2013.DATA output_dir=OPM
test $? -eq 0 || exit 1
PATH=$WORKSPACE/serial/install/bin:$PATH ./plotwells.sh

popd
