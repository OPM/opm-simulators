#!/bin/bash

OPM_DATA_ROOT=$1

# Copy results from a test run to refence dir
# $1 = source directory to copy data from
# $2 = destination directory to copy data to
# $3 = base file name for files to copy
# $4...$@ = file types to copy
copyToReferenceDir () {
  SRC_DIR=$1
  DST_DIR=$2;
  STEM=$3;
  FILETYPES=${@:4};

  for filetype in $FILETYPES; do
    cp "$SRC_DIR$STEM.$filetype" $DST_DIR
  done
}

for test_name in ${@:2}; do
  if grep -q "spe1" <<< $test_name
  then
    copyToReferenceDir \
      $configuration/build-opm-simulators/tests/results/flow_sequential+spe1/ \
      $OPM_DATA_ROOT/spe1/opm-simulation-reference/ \
      SPE1CASE1 \
      EGRID INIT SMSPEC UNRST UNSMRY

    copyToReferenceDir \
      $configuration/build-opm-simulators/tests/results/flow_sequential+spe1/ \
      $OPM_DATA_ROOT/spe1/opm-simulation-reference/ \
      SPE1CASE2 \
      EGRID INIT SMSPEC UNRST UNSMRY
  fi

  if grep -q "spe3" <<< $2
  then
    copyToReferenceDir \
      $configuration/build-opm-simulators/tests/results/flow_sequential+spe3/ \
      $OPM_DATA_ROOT/spe3/opm-simulation-reference/ \
      SPE3CASE1 \
      EGRID INIT PRT SMSPEC UNRST UNSMRY
  fi

  if grep -q "spe9" <<< $2
  then
    copyToReferenceDir \
      $configuration/build-opm-simulators/tests/results/flow+spe9/ \
      $OPM_DATA_ROOT/spe9/opm-simulation-reference/ \
      SPE9_CP_SHORT \
      EGRID INIT PRT SMSPEC UNRST UNSMRY
  fi
done





echo -e "update reference data for $2\n" > /tmp/cmsg
for dep in ert opm-common opm-core opm-grid opm-material opm-parser opm-output
do
  pushd $WORKSPACE/deps/$dep > /dev/null
  name=`printf "%-14s" $dep`
  rev=`git rev-parse HEAD`
  echo -e "$name = $rev" >> /tmp/cmsg
  popd > /dev/null
done
echo -e "opm-simulators = `git rev-parse HEAD`" >> /tmp/cmsg

cd $OPM_DATA_ROOT
git commit -a -t /tmp/cmsg
