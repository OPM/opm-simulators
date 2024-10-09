#!/bin/bash

# This runs a simulator from start to end, then a restarted run of the
# simulator. Finally the *summary output* from the the two simulations is
# compared. Primarily the driver 'run-restart-regressionTest.sh' should be used,
# but in situations where it is difficult to get the restart files to agree this
# script can be used as second-best alternative.

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
  echo -e "\t\t -d <path>     Path to restart deck tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  echo -e "\t\t -s <step>     Step to do restart testing from"
  exit 1
fi

OPTIND=1
while getopts "i:r:b:f:a:t:c:e:d:s:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    d) RST_DECK_COMMAND=${OPTARG} ;;
    s) RESTART_STEP=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

BASE_NAME=${FILENAME}_RESTART

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} --output-dir=${RESULT_PATH} ${TEST_ARGS}

test $? -eq 0 || exit 1

${RST_DECK_COMMAND} -m inline -s ${INPUT_DATA_PATH}/${FILENAME}.DATA ${FILENAME}.UNRST:${RESTART_STEP} ${BASE_NAME}.DATA

${BINPATH}/${EXE_NAME} ${BASE_NAME} --output-dir=${RESULT_PATH} ${TEST_ARGS}
test $? -eq 0 || exit 1

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_ECL_COMMAND} -R -t SMRY ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -R -t SMRY ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/${FILENAME}_RESTART ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
