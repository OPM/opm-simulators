#!/bin/bash

# This test driver just verifies that the simulator can be initialized and start simulating
# for a restarted run. Does not (yet) do any regression testing

INPUT_DATA_PATH="$1"
RESULT_PATH="$2"
BINPATH="$3"
FILENAME="$4"
SCHED_RESTART="$5"
ABS_TOL="$6"
REL_TOL="$7"
COMPARE_ECL_COMMAND="$8"
OPM_PACK_COMMAND="$9"
PARALLEL="${10}"
EXE_NAME="${11}"
shift 11
TEST_ARGS="$@"

BASE_NAME=${FILENAME}_RESTART.DATA


if test $PARALLEL -eq 1
then
  CMD_PREFIX="mpirun -np 4 "
else
  CMD_PREFIX=""
fi


rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
${OPM_PACK_COMMAND} -c ${RESULT_PATH} ${INPUT_DATA_PATH}/${FILENAME}_RESTART.DATA


cd ${RESULT_PATH}
${CMD_PREFIX} ${BINPATH}/${EXE_NAME} ${BASE_NAME} --enable-adaptive-time-stepping=false --output-dir=${RESULT_PATH} ${TEST_ARGS}
test $? -eq 0 || exit 1

ecode=0
echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l -t UNRST ${RESULT_PATH}/restart/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
    ecode=1
    ${COMPARE_ECL_COMMAND} -a -l -t UNRST ${RESULT_PATH}/restart/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
