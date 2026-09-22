#!/bin/bash

# Runs one deck in grid and aux numerical-aquifer mode and compares the summaries.
# Time steps are pinned and tolerances tightened, since DOF ordering alone changes the
# adaptive substeps and makes the runs drift apart by ~1e-3.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read deck from"
  echo -e "\t\t -f <filename> Deck file name"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -a <tol>      Absolute tolerance in comparison"
  echo -e "\t\t -t <tol>      Relative tolerance in comparison"
  echo -e "\t\t -c <path>     Path to comparison tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  echo -e "\tOptional options:"
  echo -e "\t\t -x            Compare only the summary vectors both runs produce.  Needed"
  echo -e "\t\t               for a deck that asks for block data at an aquifer cell: that"
  echo -e "\t\t               cell is a grid cell in one representation and not in the"
  echo -e "\t\t               other, so the vector exists in one run only."
  exit 1
fi

IGNORE_ONE_SIDED=""
OPTIND=1
while getopts "i:f:r:a:t:c:e:x" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    x) IGNORE_ONE_SIDED="-y" ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

PINNED_ARGS="--enable-adaptive-time-stepping=false \
             --tolerance-cnv=1e-8 \
             --tolerance-mb=1e-12 \
             --newton-min-iterations=2"

mkdir -p ${RESULT_PATH}
for MODE in grid aux
do
  rm -rf ${RESULT_PATH}/${MODE}
  mkdir -p ${RESULT_PATH}/${MODE}
  "${EXE_NAME}" ${INPUT_DATA_PATH}/${FILENAME} ${TEST_ARGS} ${PINNED_ARGS} \
                --numerical-aquifer-mode=${MODE} \
                --output-dir=${RESULT_PATH}/${MODE}
  test $? -eq 0 || exit 1
done

echo "=== Comparing the summary of the two numerical-aquifer representations ==="
if test -n "${IGNORE_ONE_SIDED}"
then
  echo "    (comparing only the vectors both runs produce -- see -x)"
fi
${COMPARE_ECL_COMMAND} -t SMRY ${IGNORE_ONE_SIDED} -a ${RESULT_PATH}/grid/${FILENAME} \
                                  ${RESULT_PATH}/aux/${FILENAME} \
                                  ${ABS_TOL} ${REL_TOL}
exit $?
