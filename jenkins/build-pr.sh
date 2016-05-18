#!/bin/bash

source `dirname $0`/build-opm-simulators.sh

declare -a upstreams
upstreams=(opm-parser
           opm-material
           opm-core
           opm-grid
           opm-output)

declare -A upstreamRev
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master
upstreamRev[opm-grid]=master
upstreamRev[opm-output]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

if grep -q "ert=" <<< $ghprbCommentBody
then
  ERT_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ert=([0-9]+).*/\1/g'`/merge
fi

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

echo "Building with ert=$ERT_REVISION opm-common=$OPM_COMMON_REVISION opm-parser=${upstreamRev[opm-parser]} opm-material=${upstreamRev[opm-material]} opm-core=${upstreamRev[opm-core]} opm-grid=${upstreamRev[opm-grid]} opm-output=${upstreamRev[opm-output]} opm-simulators=$sha1"

build_opm_simulators
test $? -eq 0 || exit 1

cp serial/build-opm-simulators/testoutput.xml .
