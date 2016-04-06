#!/bin/bash

# Clone opm-data if necessary
if ! test -d deps/opm-data
then
  OPM_DATA_REVISION="master"
  if grep -q "opm-data=" <<< $ghprbCommentBody
  then
    OPM_DATA_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-data=([0-9]+).*/\1/g'`/merge
  fi
  source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh
  clone_module opm-data $OPM_DATA_REVISION
fi

pushd .
cd deps/opm-data

# Run the norne case
cd norne
$WORKSPACE/serial/build-opm-autodiff/bin/flow deck_filename=NORNE_ATW2013.DATA output_dir=OPM
test $? -eq 0 || exit 1
PATH=$WORKSPACE/serial/install/bin:$PATH ./plotwells.sh

popd
