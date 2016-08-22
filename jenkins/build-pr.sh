#!/bin/bash

source `dirname $0`/build-opm-simulators.sh

declare -a upstreams
upstreams=(ert
           opm-parser
           opm-output
           opm-material
           opm-core
           opm-grid)

declare -A upstreamRev
upstreamRev[ert]=master
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master
upstreamRev[opm-grid]=master
upstreamRev[opm-output]=master

OPM_COMMON_REVISION=master

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  OPM_COMMON_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

for upstream in ${upstreams[*]}
do
  if grep -q "$upstream=" <<< $ghprbCommentBody
  then
    upstreamRev[$upstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*$upstream=([0-9]+).*/\1/g"`/merge
  fi
done

echo "Building with opm-common=$OPM_COMMON_REVISION ert=${upstreamRev[ert]} opm-parser=${upstreamRev[opm-parser]} opm-material=${upstreamRev[opm-material]} opm-core=${upstreamRev[opm-core]} opm-grid=${upstreamRev[opm-grid]} opm-output=${upstreamRev[opm-output]} opm-simulators=$sha1"

build_opm_simulators
test $? -eq 0 || exit 1

cp serial/build-opm-simulators/testoutput.xml .
