#!/bin/bash

# Plots the wells of the norne runs against the reference solutions.
# Includes the serial run and, if procs is given, the parallel run.
configuration=$1
procs=$2

cd $WORKSPACE/deps/opm-tests/norne

./plotwells.sh $WORKSPACE/$configuration/install/bin "ECL.2014.2 opm-simulation-reference/flow_legacy" norne-wells ${procs}
./plotwells.sh $WORKSPACE/$configuration/install/bin "opm-simulation-reference/flow_legacy" norne-wells-noecl ${procs}
