#!/bin/bash

# the script is intended to be run like this: bash hipify_file.sh ${PROJECT_BUILD_DIR} ${PROJECT_BINARY_DIR}
# it should be run automatically on the correct files through cmake
input_file=$1
output_file=$2

# make sure the output folder exists
mkdir -p $(dirname $output_file)

# hipify out-of-place
hipify-perl $input_file > $output_file

# expand includes so we only need include_directories (path to hip)
sed -i 's/^#include <hipblas\.h>/#include <hipblas\/hipblas.h>/g' $output_file
sed -i 's/^#include <hipsparse\.h>/#include <hipsparse\/hipsparse.h>/g' $output_file
# make sure includes refer to hipistl/ files (the ones that are also hipified)
sed -i 's/cuistl\//hipistl\//g' $output_file

echo "$output_file hipified"
