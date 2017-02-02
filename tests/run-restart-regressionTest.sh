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

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA timestep.adaptive=false
test $? -eq 0 || exit 1
${BINPATH}/${EXE_NAME} ${TEST_ARGS}_RESTART.DATA timestep.adaptive=false
test $? -eq 0 || exit 1

ecode=0
${COMPARE_SUMMARY_COMMAND} -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
fi

${COMPARE_ECL_COMMAND} -l ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
test $? -eq 0 || ecode=1

exit $ecode
