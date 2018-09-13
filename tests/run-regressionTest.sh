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
if test "${EXE_NAME}" = "flow"; then
    ${BINPATH}/${EXE_NAME} ${TEST_ARGS} --output-dir=${RESULT_PATH}
else
    ${BINPATH}/${EXE_NAME} ${TEST_ARGS} output_dir=${RESULT_PATH}
fi
test $? -eq 0 || exit 1
cd ..

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_ECL_COMMAND} -t SMRY ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -t SMRY ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

ignore_extra_kw=""
if grep -q "ignore_extra" <<< $ghprbCommentBody
then
    ignore_extra_kw="-x"
fi

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND}  ${ignore_extra_kw} ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} -a ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

echo "=== Executing comparison for init file ==="
${COMPARE_ECL_COMMAND} ${ignore_extra_kw} -t INIT ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} -a -t INIT ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
