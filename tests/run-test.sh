#!/bin/bash

# This simply runs a simulator.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read deck from"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -b <path>     Path to simulator binary"
  echo -e "\t\t -f <filename> Deck file name"
  echo -e "\t\t -e <filename> Simulator binary to use"
  echo -e "\t\t -n <procs >   Number of MPI processes to use"
  exit 1
fi

OPTIND=1
while getopts "i:r:b:f:e:n:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    n) MPI_PROCS=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
if (( ${MPI_PROCS} > 1))
then
  mpirun -np ${MPI_PROCS} ${BINPATH}/${EXE_NAME} ${TEST_ARGS} --output-dir=${RESULT_PATH} "${INPUT_DATA_PATH}/${FILENAME}.DATA"
else
  ${BINPATH}/${EXE_NAME} ${TEST_ARGS} --output-dir=${RESULT_PATH} "${INPUT_DATA_PATH}/${FILENAME}.DATA"
fi
test $? -eq 0 || exit 1
