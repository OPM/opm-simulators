#!/bin/bash

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
${BINPATH}/${EXE_NAME} ${TEST_ARGS}
cd ..

ecode=0
${COMPARE_SUMMARY_COMMAND} -r ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
fi

${COMPARE_ECL_COMMAND}  ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${FILENAME} ${ABS_TOL} ${REL_TOL}
test $? -eq 0 || ecode=1

${COMPARE_ECL_COMMAND} -t INIT ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${FILENAME} ${ABS_TOL} ${REL_TOL}
test $? -eq 0 || ecode=1

exit $ecode
