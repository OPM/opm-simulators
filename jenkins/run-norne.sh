#!/bin/bash

configuration=$1
procs=$2

cd $WORKSPACE/deps/opm-tests/norne

# Run the norne case
if test -n "$procs"
then
  mpirun -np $procs $WORKSPACE/$configuration/build-opm-simulators/bin/flow --output-dir=flow_${procs}_proc NORNE_ATW2013.DATA
else
  $WORKSPACE/$configuration/build-opm-simulators/bin/flow --output-dir=flow NORNE_ATW2013.DATA
fi
test $? -eq 0 || exit 1
./plotwells.sh $WORKSPACE/$configuration/install/bin "ECL.2014.2 opm-simulation-reference/flow_legacy" norne-wells
./plotwells.sh $WORKSPACE/$configuration/install/bin "opm-simulation-reference/flow_legacy" norne-wells-noecl
