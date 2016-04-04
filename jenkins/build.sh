#!/bin/bash

source `dirname $0`/build-opm-material.sh

ERT_REVISION=master
OPM_COMMON_REVISION=master
OPM_PARSER_REVISION=master

build_opm_material
test $? -eq 0 || exit 1

cp serial/build-opm-material/testoutput.xml .
