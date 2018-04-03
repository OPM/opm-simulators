#!/bin/bash

pushd .
cd deps/opm-tests
test -z $SIM && SIM=flow

# Run the norne case
cd norne
mkdir $SIM
$WORKSPACE/$configuration/build-opm-simulators/bin/$SIM deck_filename=NORNE_ATW2013.DATA output_dir=$SIM
test $? -eq 0 || exit 1
./plotwells.sh $WORKSPACE/$configuration/install/bin "ECL.2014.2 opm-simulation-reference/flow_legacy" norne-wells
./plotwells.sh $WORKSPACE/$configuration/install/bin "opm-simulation-reference/flow_legacy" norne-wells-noecl

popd
