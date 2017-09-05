#!/bin/bash

#Print commands as they execute
#set -x

declare -a upstreams
upstreams=(opm-common
           libecl
           opm-parser
           opm-output
           opm-material
           opm-grid
           opm-core
           ewoms)

declare -A upstreamRev
upstreamRev[opm-common]=master
upstreamRev[libecl]=master
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master
upstreamRev[opm-grid]=master
upstreamRev[opm-output]=master
upstreamRev[ewoms]=master

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  upstreamRev[opm-common]=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

# No downstreams currently
declare -a downstreams
declare -A downstreamRev

# Clone opm-common
pushd .
if [ ! -d $WORKSPACE/deps/opm-common ];then
    mkdir -p $WORKSPACE/deps/opm-common
    cd $WORKSPACE/deps/opm-common
    git init .
    git remote add origin https://github.com/OPM/opm-common
    git fetch --depth 1 origin ${upstreamRev[opm-common]}:branch_to_build
    test $? -eq 0 || exit 1
    git checkout branch_to_build
else
    cd $WORKSPACE/deps/opm-common
    #git fetch --depth 1 origin ${upstreamRev[opm-common]}:branch_to_build
    #test $? -eq 0 || exit 1
    git checkout branch_to_build   
    git pull
fi
       
popd

source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

parseRevisions
printHeader opm-simulators

# Setup opm-data
source $WORKSPACE/deps/opm-common/jenkins/setup-opm-data.sh

build_module_full opm-simulators
