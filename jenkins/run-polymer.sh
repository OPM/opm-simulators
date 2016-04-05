#!/bin/bash

# Build flow_polymer binary
pushd .
cd serial/build-opm-autodiff
cmake --build . --target flow_polymer
popd

# Clone opm-data if necessary
pushd .
cd deps
if ! test -d opm-data
then
  git clone --depth 1 --single-branch -b master https://github.com/OPM/opm-data
fi
cd opm-data

# Run the simple2D polymer case
cd polymer_test_suite/simple2D
$WORKSPACE/serial/build-opm-autodiff/bin/flow_polymer run.param
test $? -eq 0 || exit 1
cd ../..

# Compare OPM with eclipse reference
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py polymer_test_suite/simple2D/eclipse-simulation/ polymer_test_suite/simple2D/opm-simulation/ 2D_THREEPHASE_POLY_HETER 0.0006 0.004
test $? -eq 0 || exit 1

popd
