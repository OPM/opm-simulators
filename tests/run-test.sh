#!/bin/bash

# This simply runs a simulator.

INPUT_DATA_PATH="$1"
RESULT_PATH="$2"
BINPATH="$3"
EXE_NAME="$4"
FILENAME="$5"
MPI_PROCS="$6"
shift 7
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
if (( ${MPI_PROCS} > 1))
then
  mpirun -np ${MPI_PROCS} ${BINPATH}/${EXE_NAME} ${TEST_ARGS} --output-dir=${RESULT_PATH}
else
  ${BINPATH}/${EXE_NAME} ${TEST_ARGS} --output-dir=${RESULT_PATH}
fi
test $? -eq 0 || exit 1
