#!/bin/bash

# Runs a deck twice with the same simulator, once with the intensive
# quantities computed on the CPU and once on the GPU
# (--experimental-compute-properties-on-gpu=true), then compares the summary
# and restart files of the two runs.
# Meant to track regressions in the GPU intensive quantities.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read deck from"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -f <filename> Deck file name"
  echo -e "\t\t -a <tol>      Absolute tolerance in comparison"
  echo -e "\t\t -t <tol>      Relative tolerance in comparison"
  echo -e "\t\t -c <path>     Path to comparison tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  exit 1
fi

OPTIND=1
while getopts "i:r:f:a:t:c:e:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
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
mkdir -p ${RESULT_PATH}/cpu ${RESULT_PATH}/gpu

"${EXE_NAME}" ${INPUT_DATA_PATH}/${FILENAME} ${TEST_ARGS} \
  --experimental-compute-properties-on-gpu=false --output-dir=${RESULT_PATH}/cpu
test $? -eq 0 || exit 1

"${EXE_NAME}" ${INPUT_DATA_PATH}/${FILENAME} ${TEST_ARGS} \
  --experimental-compute-properties-on-gpu=true --output-dir=${RESULT_PATH}/gpu
test $? -eq 0 || exit 1

# Without this line the GPU path could silently fall back to the CPU and
# the comparison would pass trivially.
grep -q "GpuBlackoilIntensiveQuantitiesDispatcher\] initialized" ${RESULT_PATH}/gpu/${FILENAME}.PRT \
  || { echo "The GPU run did not use the GPU intensive-quantities dispatcher"; exit 1; }

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_ECL_COMMAND} -t SMRY -R ${RESULT_PATH}/cpu/${FILENAME} ${RESULT_PATH}/gpu/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -t SMRY -a -R ${RESULT_PATH}/cpu/${FILENAME} ${RESULT_PATH}/gpu/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l -t UNRST ${RESULT_PATH}/cpu/${FILENAME} ${RESULT_PATH}/gpu/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} -a -l -t UNRST ${RESULT_PATH}/cpu/${FILENAME} ${RESULT_PATH}/gpu/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
