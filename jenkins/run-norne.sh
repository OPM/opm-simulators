#!/bin/bash

pushd .
cd deps/opm-tests
test -z $SIM && SIM=flow

# Run the norne case
cd norne
if test -n "$1"
then
  mpirun -np $1 $WORKSPACE/$configuration/build-opm-simulators/bin/$SIM --output-dir=${SIM}_${1}_proc NORNE_ATW2013.DATA
else
  $WORKSPACE/$configuration/build-opm-simulators/bin/$SIM --output-dir=$SIM NORNE_ATW2013.DATA
fi
test $? -eq 0 || exit 1
./plotwells.sh $WORKSPACE/$configuration/install/bin "ECL.2014.2 opm-simulation-reference/flow_legacy" norne-wells
./plotwells.sh $WORKSPACE/$configuration/install/bin "opm-simulation-reference/flow_legacy" norne-wells-noecl

popd
