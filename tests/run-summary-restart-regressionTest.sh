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
  echo -e "\t\t -p <path>     Path to deck packing tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  exit 1
fi

OPTIND=1
while getopts "i:r:b:f:a:t:c:p:e:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    p) OPM_PACK_COMMAND=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

BASE_NAME=${FILENAME}_RESTART.DATA

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} --output-dir=${RESULT_PATH} ${TEST_ARGS}

test $? -eq 0 || exit 1

${OPM_PACK_COMMAND} -o ${BASE_NAME} ${INPUT_DATA_PATH}/${FILENAME}_RESTART.DATA

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
