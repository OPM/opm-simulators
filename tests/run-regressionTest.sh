#!/bin/bash

# This runs a simulator, then compares the summary, restart and init
# files against a reference.

INPUT_DATA_PATH="$1"
RESULT_PATH="$2"
BINPATH="$3"
FILENAME="$4"
ABS_TOL="$5"
REL_TOL="$6"
COMPARE_ECL_COMMAND="$7"
EXE_NAME="${8}"
shift 8
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1
cd ..


ecode=0
echo "=== Executing comparison for EGRID, INIT, UNRST and RFT files if these exists in refrece folder ==="
${COMPARE_ECL_COMMAND}  ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${RESULT_PATH}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a  ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${RESULT_PATH}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
