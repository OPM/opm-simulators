#!/bin/bash

source `dirname $0`/build-opm-material.sh

declare -a upstreams
upstreams=(opm-parser)

declare -A upstreamRev
upstreamRev[opm-parser]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

build_opm_material
test $? -eq 0 || exit 1

cp serial/build-opm-material/testoutput.xml .
