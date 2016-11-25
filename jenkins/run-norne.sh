#!/bin/bash

pushd .
cd deps/opm-data
test -z $SIM && SIM=flow

# Run the norne case
cd norne
$WORKSPACE/$configuration/build-opm-simulators/bin/$SIM deck_filename=NORNE_ATW2013.DATA output_dir=OPM
test $? -eq 0 || exit 1
LD_LIBRARY_PATH=$WORKSPACE/$configuration/lib/x86_64-linux-gnu PATH=$WORKSPACE/$configuration/install/bin:$PATH ./plotwells.sh

popd
