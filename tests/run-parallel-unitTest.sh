#!/bin/bash
# This executes a unit test in parallel.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -n <procs>    Number of MPI Processes to use"
  echo -e "\t\t -e <filename> Simulator binary to use"
  exit 1
fi

OPTIND=1
while getopts "n:e:" OPT
do
  case "${OPT}" in
    e) EXE_NAME=${OPTARG} ;;
    n) NP=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

mpirun -np $NP "${EXE_NAME}" ${TEST_ARGS}
