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
COMPARE_SUMMARY_COMMAND="$7"
COMPARE_ECL_COMMAND="$8"
OPM_PACK_COMMAND="$9"
PARALLEL="${10}"
EXE_NAME="${11}"
shift 11
TEST_ARGS="$@"

BASE_NAME=`basename ${TEST_ARGS}_RESTART.DATA`

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
if test $PARALLEL -eq 1
then
  CMD_PREFIX="mpirun -np 4 "
else
  CMD_PREFIX=""
fi
if test "${EXE_NAME}" = "flow"; then
    ${CMD_PREFIX} ${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA --enable-adaptive-time-stepping=false --output-dir=${RESULT_PATH}
else
    ${CMD_PREFIX} ${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA timestep.adaptive=false output_dir=${RESULT_PATH}
fi

test $? -eq 0 || exit 1

${OPM_PACK_COMMAND} -o ${BASE_NAME} ${TEST_ARGS}_RESTART.DATA

if test "${EXE_NAME}" = "flow"; then
    ${CMD_PREFIX} ${BINPATH}/${EXE_NAME} ${BASE_NAME} --enable-adaptive-time-stepping=false --output-dir=${RESULT_PATH}
else
    ${CMD_PREFIX} ${BINPATH}/${EXE_NAME} ${BASE_NAME} timestep.adaptive=false output_dir=${RESULT_PATH}
fi
test $? -eq 0 || exit 1

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_SUMMARY_COMMAND} -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_SUMMARY_COMMAND} -a -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
fi

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -l ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
