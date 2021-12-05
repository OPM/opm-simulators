#!/bin/bash

# This runs a simulator, then compares the summary, restart and init
# files against a reference.

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
  echo -e "\tOptional options:"
  echo -e "\t\t -s <step>     Step to do restart testing from"
  echo -e "\t\t -h value      sched_restart value to use in restart test"
  exit 1
fi

RESTART_STEP=""
OPTIND=1
while getopts "i:r:b:f:a:t:c:d:s:e:h:" OPT
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
    h) RESTART_SCHED=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1
cd ..


ecode=0

ignore_extra_kw=""
if grep -q "ignore_extra" <<< $ghprbCommentBody
then
    ignore_extra_kw="-x"
fi

echo "=== Executing comparison for EGRID, INIT, UNRST and RFT files if these exists in reference folder ==="
${COMPARE_ECL_COMMAND} ${ignore_extra_kw} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${RESULT_PATH}/${FILENAME} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
  ecode=1
  ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} -a  ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${FILENAME} ${RESULT_PATH}/${FILENAME} ${ABS_TOL} ${REL_TOL}
fi

RSTEPS=(${RESTART_STEP//,/ })

for STEP in "${RSTEPS[@]}"
do
  echo "=== Executing restart run from step: ${STEP} ==="
  mkdir -p ${RESULT_PATH}/restart
  cp -f ${RESULT_PATH}/${FILENAME}.UNRST ${RESULT_PATH}/restart
  ${RST_DECK_COMMAND}  ${INPUT_DATA_PATH}/${FILENAME}.DATA ${FILENAME}:${STEP} -m inline -s > ${RESULT_PATH}/restart/${FILENAME}_RESTART_${STEP}.DATA
  cd ${RESULT_PATH}/restart
  if test -n "$RESTART_SCHED"
  then
    sched_rst="--sched-restart=${RESTART_SCHED}"
  fi
  ${BINPATH}/${EXE_NAME} ${TEST_ARGS} ${sched_rst} --output-dir=${RESULT_PATH}/restart ${FILENAME}_RESTART_${STEP}
  test $? -eq 0 || exit 1

  echo "=== Executing comparison for EGRID, INIT, UNRST and RFT files for restarted run ==="
  ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/restart/${FILENAME}_RESTART_${STEP} ${RESULT_PATH}/restart/${FILENAME}_RESTART_${STEP} ${ABS_TOL} ${REL_TOL}
  if [ $? -ne 0 ]
  then
    ecode=1
    ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} -a  ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/restart/${FILENAME}_RESTART_${STEP} ${RESULT_PATH}/restart/${FILENAME}_RESTART_${STEP} ${ABS_TOL} ${REL_TOL}
  fi
done

exit $ecode
