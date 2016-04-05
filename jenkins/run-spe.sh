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

# Run the SPE1/3/9 cases
cd spe1
$WORKSPACE/serial/build-opm-autodiff/bin/flow deck_filename=SPE1CASE2.DATA
test $? -eq 0 || exit 1
cd ..
cd spe3
$WORKSPACE/serial/build-opm-autodiff/bin/flow max_iter=50 deck_filename=SPE3CASE1.DATA
test $? -eq 0 || exit 1
cd ..
cd spe9
$WORKSPACE/serial/build-opm-autodiff/bin/flow max_iter=50 deck_filename=SPE9_CP.DATA
test $? -eq 0 || exit 1
cd ..

# Compare OPM with eclipse reference
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py spe1/eclipse-simulation/ spe1/ SPE1CASE2 0.01 0.01
test $? -eq 0 || exit 1
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py spe3/eclipse-simulation/ spe3/ SPE3CASE1 0.02 0.02
test $? -eq 0 || exit 1
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py spe9/eclipse-simulation/ spe9/ SPE9_CP 0.002 0.001
test $? -eq 0 || exit 1

# Compare OPM with OPM reference
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py spe1/opm-simulation-reference/ spe1/ SPE1CASE2 0.001 0.001
test $? -eq 0 || exit 1
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py spe3/opm-simulation-reference/ spe3/ SPE3CASE1 0.001 0.001
test $? -eq 0 || exit 1
PYTHONPATH=$WORKSPACE/serial/install/lib/python2.7/dist-packages/ python output_comparator/src/compare_eclipse.py spe9/opm-simulation-reference/ spe9/ SPE9_CP 0.002 0.007
test $? -eq 0 || exit 1

popd
