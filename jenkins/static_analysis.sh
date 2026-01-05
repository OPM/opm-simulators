#!/bin/bash

declare -a upstreams
upstreams=(opm-common
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-common]=master
upstreamRev[opm-grid]=master

if grep -q "opm-common=" <<< $ghprbCommentBody
then
    if test -n "$absolute_revisions"
    then
        upstreamRev[opm-common]=$(echo $ghprbCommentBody | sed -r 's/.*opm-common=([^ ]+)/\1/g')
    else
        upstreamRev[opm-common]=pull/$(echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g')/merge
    fi
fi

# Clone opm-common
mkdir -p $WORKSPACE/deps/opm-common
pushd $WORKSPACE/deps/opm-common
repo_root=${OPM_REPO_ROOT:-https://github.com/OPM}
git init .
git remote add origin ${repo_root}/opm-common
git fetch --depth 1 origin ${upstreamRev[opm-common]}:branch_to_build
test $? -eq 0 || exit 1
git checkout branch_to_build
popd

source ${TOOLCHAIN_DIR}/build-configurations-sca.sh
source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

$WORKSPACE/jenkins/build.sh

run_static_analysis opm-simulators
