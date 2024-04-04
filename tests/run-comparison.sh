#!/bin/bash

# This does two simulator runs and compares the summary files against each other.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read first deck from"
  echo -e "\t\t -j <path>     Path to read second deck from"
  echo -e "\t\t -f <filename> First deck file name"
  echo -e "\t\t -g value      Second deck file name"
  echo -e "\t\t -y value      Ignore extra keywords in the run. For -y 'BOTH', extra keywords in both runs will be ignored. For -y 'SECOND', extra keywords in the second deck will be ignored."
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -b <path>     Path to simulator binary"
  echo -e "\t\t -a <tol>      Absolute tolerance in comparison"
  echo -e "\t\t -t <tol>      Relative tolerance in comparison"
  echo -e "\t\t -c <path>     Path to comparison tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  exit 1
fi

RESTART_STEP=""
OPTIND=1
while getopts "i:j:f:g:r:b:a:t:c:e:y:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH1=${OPTARG} ;;
    j) INPUT_DATA_PATH2=${OPTARG} ;;
    f) FILENAME1=${OPTARG} ;;
    g) FILENAME2=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    y) IGNORE_EXTRA_KW=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH1}/${FILENAME1} ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1
${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH2}/${FILENAME2} ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1
cd ..


ecode=0

ignore_extra_kw=""
if grep -q "SECOND" <<< $IGNORE_EXTRA_KW
then
    ignore_extra_kw="-x"
fi

if grep -q "BOTH" <<< $IGNORE_EXTRA_KW
then
  ignore_extra_kw="-y"
fi

echo "=== Executing comparison for EGRID, INIT, UNRST and RFT files if these exists in reference folder ==="
${COMPARE_ECL_COMMAND} -t SMRY ${ignore_extra_kw} ${RESULT_PATH}/${FILENAME1} ${RESULT_PATH}/${FILENAME2} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -t SMRY ${ignore_extra_kw} -a ${RESULT_PATH}/${FILENAME1} ${RESULT_PATH}/${FILENAME2} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
