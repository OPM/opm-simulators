#!/bin/bash

# This runs two parallel cases, one of them with one more process
# and the --test-split-communicator=true option,
# then compares the summary and restart files from the two runs.
# Meant to track regression of the treatment of MPI communicators.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read deck from"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -b <path>     Path to simulator binary"
  echo -e "\t\t -f <filename> Deck file name"
  echo -e "\t\t -a <tol>      Absolute tolerance in comparison"
  echo -e "\t\t -t <tol>      Relative tolerance in comparison"
  echo -e "\t\t -c <path>     Path to comparison tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  exit 1
fi

BASE_MPI_PROCS=3
TEST_MPI_PROCS=4 # should be 1 more than the base
OPTIND=1
while getopts "i:r:b:f:a:t:c:e:n:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}

echo mpirun -np ${BASE_MPI_PROCS} ${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME}.DATA ${TEST_ARGS} --output-dir=${RESULT_PATH}/base
mpirun -np ${BASE_MPI_PROCS} ${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME}.DATA ${TEST_ARGS} --output-dir=${RESULT_PATH}/base
test $? -eq 0 || exit 1

echo mpirun -np ${TEST_MPI_PROCS} ${BINPATH}/${EXE_NAME} --test-split-communicator=true ${INPUT_DATA_PATH}/${FILENAME}.DATA ${TEST_ARGS} --output-dir=${RESULT_PATH}/test
mpirun -np ${TEST_MPI_PROCS} ${BINPATH}/${EXE_NAME} --test-split-communicator=true ${INPUT_DATA_PATH}/${FILENAME}.DATA ${TEST_ARGS} --output-dir=${RESULT_PATH}/test
test $? -eq 0 || exit 1

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_ECL_COMMAND} -t SMRY -R ${RESULT_PATH}/base/${FILENAME} ${RESULT_PATH}/test/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -t SMRY -a -R ${RESULT_PATH}/base/${FILENAME} ${RESULT_PATH}/test/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l -t UNRST ${RESULT_PATH}/base/${FILENAME} ${RESULT_PATH}/test/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -l -t UNRST ${RESULT_PATH}/base/${FILENAME} ${RESULT_PATH}/test/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
