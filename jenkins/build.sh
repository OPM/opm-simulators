#!/bin/bash

source `dirname $0`/build-opm-material.sh

declare -a upstreams
upstreams=(ert
           opm-parser)

declare -A upstreamRev
upstreamRev[ert]=master
upstreamRev[opm-parser]=master

OPM_COMMON_REVISION=master

build_opm_material
test $? -eq 0 || exit 1

cp serial/build-opm-material/testoutput.xml .
