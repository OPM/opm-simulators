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

build_opm_simulators
test $? -eq 0 || exit 1

cp serial/build-opm-simulators/testoutput.xml .
