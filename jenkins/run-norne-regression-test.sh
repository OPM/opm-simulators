#!/bin/bash 
set -e 



SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
DIR=$(dirname "$SCRIPT")
cd $DIR

# Downloads some of the opm repos from github. 
# Installs them as INSTALL_PREFIX
#  
./build.sh 

OPM_DATA_PATH=${DIR}/../../opm-data/norne

RESULT_PATH=${DIR}/../build/tests

BINPATH=${DIR}/../build/bin
EXE_NAME=flow
TEST_ARG=${OPM_DATA_PATH}/NORNE_ATW2013.DATA

rm -Rf  ${RESULT_PATH};
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${TEST_ARG}
cd ..

${DIR}/../../opm-output/build/bin/compareSummary -r ${RESULT_PATH}/NORNE_ATW2013 ${OPM_DATA_PATH}/OPM/opm-solution-reference/NORNE_ATW2013 1e-6 1e-4
${DIR}/../../opm-output/build/bin/compareECL -t INIT ${RESULT_PATH}/NORNE_ATW2013 ${OPM_DATA_PATH}/OPM/opm-solution-reference/NORNE_ATW2013 1e-6 1e-4