#!/bin/bash

pushd .
cd deps/opm-data

# Run the simple2D polymer case
cd polymer_simple2D
$WORKSPACE/$configuration/build-opm-simulators/bin/flow_polymer run.param
test $? -eq 0 || exit 1
cd ..

# Compare OPM with eclipse reference
compareECL=$WORKSPACE/$configuration/install/bin/compareECL
reffile=polymer_simple2D/eclipse-simulation/2D_THREEPHASE_POLY_HETER
opmfile=polymer_simple2D/opm-simulation-reference/flow_polymer/2D_THREEPHASE_POLY_HETER
$compareECL $reffile $opmfile 1.0 0.004 -k SGAS
test $? -eq 0 || exit 1
$compareECL $reffile $opmfile 1.0 0.004 -k SWAT
test $? -eq 0 || exit 1
$compareECL $reffile $opmfile 0.0006 1.0 -k PRESSURE
test $? -eq 0 || exit 1

popd
