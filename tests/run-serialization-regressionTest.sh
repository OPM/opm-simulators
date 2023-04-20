#!/bin/bash

# This runs a simulator from start to end, then a restarted
# run of the simulator, before comparing the output from the two runs.
# This is meant to track regressions in the restart support.

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
  echo -e "\t\t -s <step>     Step to do restart testing from"
  exit 1
fi

OPTIND=1
MPI_PROCS=1
while getopts "i:r:b:f:a:t:c:e:n:d:s:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    d) ;;
    s) RESTART_STEP=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    n) MPI_PROCS=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

BASE_NAME=${FILENAME}
if test $MPI_PROCS -gt 1
then
  CMD_PREFIX="mpirun -np $MPI_PROCS "
fi

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${CMD_PREFIX}${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} --output-dir=${RESULT_PATH} ${TEST_ARGS} --save-step=${RESTART_STEP}

test $? -eq 0 || exit 1

mkdir -p ${RESULT_PATH}/restart
${CMD_PREFIX}${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} --output-dir=${RESULT_PATH}/restart ${TEST_ARGS} --load-step=${RESTART_STEP} --save-file=${RESULT_PATH}/${FILENAME}.OPMRST
test $? -eq 0 || exit 1

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l -t UNRST ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/restart/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -l -t UNRST ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/restart/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
