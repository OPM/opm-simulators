#!/bin/bash

# This performs a serial and a parallel for a simulator,
# then compares the summary and restart files from the two runs.
# Meant to track regression in parallel simulators.

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

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
if test "${EXE_NAME}" = "flow"; then
    ${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --output-dir=${RESULT_PATH}
else
    ${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA linear_solver_reduction=1e-7 tolerance_cnv=5e-6 tolerance_mb=1e-8 output_dir=${RESULT_PATH}
fi

test $? -eq 0 || exit 1
mkdir mpi
cd mpi
if test "${EXE_NAME}" = "flow"; then
    mpirun -np 4 ${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --output-dir=${RESULT_PATH}/mpi
else
    mpirun -np 4 ${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA linear_solver_reduction=1e-7 tolerance_cnv=5e-6 tolerance_mb=1e-8 output_dir=${RESULT_PATH}/mpi
fi
test $? -eq 0 || exit 1
cd ..

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_ECL_COMMAND} -t SMRY -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -t SMRY -a -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -l ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
