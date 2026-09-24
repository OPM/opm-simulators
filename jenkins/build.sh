#!/bin/bash

#Print commands as they execute
#set -x

declare -a upstreams
upstreams=(opm-common
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-common]=master
upstreamRev[opm-grid]=master

# No downstreams currently
declare -a downstreams
declare -A downstreamRev

# Fetch opm-common before loading its shared Jenkins helpers.
source "$WORKSPACE/jenkins/checkout-opm-common.sh"

source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

parseRevisions
printHeader opm-simulators

clone_repositories opm-simulators

# Setup opm-data
source $WORKSPACE/deps/opm-common/jenkins/setup-opm-tests.sh

build_module_full opm-simulators
