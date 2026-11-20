#!/bin/bash

# Runs a single porosity deck and a dual-continuum deck, and compares the cells of one
# continuum of the dual-continuum run with the single porosity run.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read the single porosity deck from"
  echo -e "\t\t -j <path>     Path to read the dual-continuum deck from"
  echo -e "\t\t -f <filename> Single porosity deck file name"
  echo -e "\t\t -g <filename> Dual-continuum deck file name"
  echo -e "\t\t -k <name>     Continuum to compare, fracture or matrix"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -p <path>     Path to the comparison script"
  echo -e "\t\t -e <filename> Simulator binary to use"
  exit 1
fi

OPTIND=1
while getopts "i:j:f:g:k:r:p:e:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH1=${OPTARG} ;;
    j) INPUT_DATA_PATH2=${OPTARG} ;;
    f) FILENAME1=${OPTARG} ;;
    g) FILENAME2=${OPTARG} ;;
    k) CONTINUUM=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    p) COMPARISON_SCRIPT=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
"${EXE_NAME}" ${INPUT_DATA_PATH1}/${FILENAME1} ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1
"${EXE_NAME}" ${INPUT_DATA_PATH2}/${FILENAME2} ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1

python3 -I "${COMPARISON_SCRIPT}" ${RESULT_PATH} ${FILENAME1%.DATA} ${RESULT_PATH} ${FILENAME2%.DATA} ${CONTINUUM}
