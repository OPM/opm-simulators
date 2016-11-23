#!/bin/bash
set -e

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
${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA linear_solver_reduction=1e-7 tolerance_cnv=5e-6 tolerance_mb=1e-8
mkdir mpi
cd mpi
mpirun -np 4 ${BINPATH}/${EXE_NAME} ${TEST_ARGS}.DATA linear_solver_reduction=1e-7 tolerance_cnv=5e-6 tolerance_mb=1e-8
cd ..

${COMPARE_SUMMARY_COMMAND} -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
${COMPARE_ECL_COMMAND} -l ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
