#!/bin/bash

source `dirname $0`/build-opm-autodiff.sh

ERT_REVISION=master
OPM_COMMON_REVISION=master
OPM_PARSER_REVISION=master
OPM_MATERIAL_REVISION=master
OPM_CORE_REVISION=master
OPM_GRID_REVISION=master
OPM_OUTPUT_REVISION=master

build_opm_autodiff
test $? -eq 0 || exit 1

cp serial/build-opm-autodiff/testoutput.xml .
