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
RST_DECK_COMMAND="$8"
RESTART_STEP="${9}"
RESTART_SCHED="${10}"
EXE_NAME="${11}"
shift 11
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1
cd ..


ecode=0

ignore_extra_kw=""
if grep -q "ignore_extra" <<< $ghprbCommentBody
then
    ignore_extra_kw="-x"
fi

echo "=== Executing comparison for EGRID, INIT, UNRST and RFT files if these exists in reference folder ==="
${COMPARE_ECL_COMMAND} ${ignore_extra_kw} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${RESULT_PATH}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} -a  ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${RESULT_PATH}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

if test $RESTART_STEP -ne 0
then
  echo "=== Executing restart run ==="
  if [ "$RESTART_SCHED" = "--" ]; then
      sched_rst=""
  else
      sched_rst="${RESTART_SCHED}"
  fi

  mkdir -p ${RESULT_PATH}/restart
  cp -f ${RESULT_PATH}/${FILENAME}.UNRST ${RESULT_PATH}/restart
  ${RST_DECK_COMMAND}  ${INPUT_DATA_PATH}/${FILENAME}.DATA ${FILENAME}:${RESTART_STEP} -m inline -s > ${RESULT_PATH}/restart/${FILENAME}.DATA
  cd ${RESULT_PATH}/restart
  echo ${BINPATH}/${EXE_NAME} ${TEST_ARGS} ${sched_rst} --output-dir=${RESULT_PATH}/restart ${FILENAME}
  ${BINPATH}/${EXE_NAME} ${TEST_ARGS} ${sched_rst} --output-dir=${RESULT_PATH}/restart ${FILENAME}
  test $? -eq 0 || exit 1

  echo "=== Executing comparison for EGRID, INIT, UNRST and RFT files for restarted run ==="
  ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/restart/${FILENAME} ${RESULT_PATH}/restart/${FILENAME} ${ABS_TOL} ${REL_TOL}
  if [ $? -ne 0 ]
  then
    ecode=1
    ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} -a  ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/restart/${FILENAME} ${RESULT_PATH}/restart/${FILENAME} ${ABS_TOL} ${REL_TOL}
  fi
fi

exit $ecode
