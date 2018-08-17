#!/bin/bash

function dotest {
  $compareECL $reffile $opmfile 1.0 $1 -k SGAS
  test $? -eq 0 || exit 1
  $compareECL $reffile $opmfile $2 1.0 -k SWAT
  test $? -eq 0 || exit 1
  $compareECL $reffile $opmfile $2 1.0 -k PRESSURE
  test $? -eq 0 || exit 1
}

pushd .
cd deps/opm-tests

EXE=flow_legacy

# Run the SPE1/3/9 cases
cd spe1
$WORKSPACE/$configuration/build-opm-simulators/bin/${EXE} SPE1CASE2.DATA
test $? -eq 0 || exit 1
cd ..
cd spe3
$WORKSPACE/$configuration/build-opm-simulators/bin/${EXE} --flow-newton-max-iterations=50 SPE3CASE1.DATA
test $? -eq 0 || exit 1
cd ..
cd spe9
$WORKSPACE/$configuration/build-opm-simulators/bin/${EXE} --flow-newton-max-iterations=50 SPE9_CP.DATA
test $? -eq 0 || exit 1
cd ..

compareECL=$WORKSPACE/$configuration/install/bin/compareECL

# Compare OPM with eclipse reference
reffile=spe1/eclipse-simulation/SPE1CASE2
opmfile=spe1/SPE1CASE2
dotest 0.01 0.01
test $? -eq 0 || exit 1

reffile=spe3/eclipse-simulation/SPE3CASE1
opmfile=spe3/SPE3CASE1
dotest 0.02 0.02
test $? -eq 0 || exit 1

reffile=spe9/eclipse-simulation/SPE9_CP
opmfile=spe9/SPE9_CP
dotest 0.002 0.001
test $? -eq 0 || exit 1

# Compare OPM with OPM reference
reffile=spe1/opm-simulation-reference/${EXE}/SPE1CASE2
opmfile=spe1/SPE1CASE2
dotest 0.001 0.001
test $? -eq 0 || exit 1

reffile=spe3/opm-simulation-reference/${EXE}/SPE3CASE1
opmfile=spe3/SPE3CASE1
dotest 0.001 0.001
test $? -eq 0 || exit 1

reffile=spe9/opm-simulation-reference/${EXE}/SPE9_CP
opmfile=spe9/SPE9_CP
dotest 0.002 0.007
test $? -eq 0 || exit 1

popd
