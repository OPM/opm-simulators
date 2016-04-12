#!/bin/bash

pushd .
cd deps/opm-data

# Run the simple2D polymer case
cd polymer_test_suite/simple2D
$WORKSPACE/serial/build-opm-simulators/bin/flow_polymer run.param
test $? -eq 0 || exit 1
cd ../..

# Compare OPM with eclipse reference
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py polymer_test_suite/simple2D/eclipse-simulation/ polymer_test_suite/simple2D/opm-simulation/ 2D_THREEPHASE_POLY_HETER 0.0006 0.004
test $? -eq 0 || exit 1

popd
