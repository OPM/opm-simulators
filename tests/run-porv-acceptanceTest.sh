#!/bin/bash

# This runs the initialization step of a simulator,
# then compares the resulting INIT file against a reference.
# This is meant to track regressions in INIT file writing.
# Useful for models that are too large to do simulation on
# as a regression test.

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

OPTIND=1
while getopts "i:r:b:f:a:t:c:e:d:h:" OPT
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
    d) ;; # Ignored
    h) ;; # Ignored
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

rm -Rf  ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${TEST_ARGS} --enable-dry-run=true --output-dir=${RESULT_PATH} ${INPUT_DATA_PATH}/${FILENAME}
cd ..

ecode=0
${COMPARE_ECL_COMMAND} -t INIT -k PORV ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-porevolume-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -t INIT -k PORV ${RESULT_PATH}/${FILENAME} ${INPUT_DATA_PATH}/opm-porevolume-reference/${EXE_NAME}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
