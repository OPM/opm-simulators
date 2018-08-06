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
COMPARE_SUMMARY_COMMAND="$7"
COMPARE_ECL_COMMAND="$8"
EXE_NAME="${9}"
shift 9
TEST_ARGS="$@"

rm -Rf  ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
if test "${EXE_NAME}" = "flow"; then
    ${BINPATH}/${EXE_NAME} ${TEST_ARGS} --enable-dry-run=true --ecl-output-dir=${RESULT_PATH}
else
    ${BINPATH}/${EXE_NAME} ${TEST_ARGS} nosim=true output_dir=${RESULT_PATH}
fi
cd ..

ecode=0
${COMPARE_ECL_COMMAND} -t INIT ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -t INIT ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
