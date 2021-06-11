#!/bin/bash

# This runs a simulator from start to end, then a restarted
# run of the simulator, before comparing the output from the two runs.
# This is meant to track regressions in the restart support.

INPUT_DATA_PATH="$1"
RESULT_PATH="$2"
BINPATH="$3"
FILENAME="$4"
ABS_TOL="$5"
REL_TOL="$6"
COMPARE_ECL_COMMAND="$7"
OPM_PACK_COMMAND="$8"
EXE_NAME="${9}"
shift 9
TEST_ARGS="$@"

BASE_NAME=${FILENAME}_RESTART.DATA

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} --output-dir=${RESULT_PATH} ${TEST_ARGS}

test $? -eq 0 || exit 1

${OPM_PACK_COMMAND} -o ${BASE_NAME} ${INPUT_DATA_PATH}/${FILENAME}_RESTART.DATA

${BINPATH}/${EXE_NAME} ${BASE_NAME} --output-dir=${RESULT_PATH} ${TEST_ARGS}
test $? -eq 0 || exit 1

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_ECL_COMMAND} -R -t SMRY ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -R -t SMRY ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
fi

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l -t UNRST ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -l -t UNRST ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
