#!/bin/bash

# This runs the initialization step of a simulator,
# then compares the resulting INIT file against a reference.
# This is meant to track regressions in INIT file writing.
# Useful for models that are too large to do simulation on
# as a regression test.

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

rm -Rf  ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${TEST_ARGS} --enable-dry-run=true --output-dir=${RESULT_PATH}
cd ..


ecode=0
${COMPARE_ECL_COMMAND} -t INIT -k PORV ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-porevolume-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -t INIT -k PORV ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-porevolume-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
