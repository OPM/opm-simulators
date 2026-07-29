#!/bin/bash

# Runs one deck twice -- once with each representation of its numerical aquifers -- and
# compares the two against each other.
#
# The aquifer is the same discrete system either way: the same unknowns, pore volumes,
# depths, regions and connection transmissibilities.  Only where the unknown lives
# changes, so the two runs should agree to the level of arithmetic-ordering roundoff, and
# anything looser is a real difference that wants explaining rather than a wider band.
#
# That comparison is only meaningful with the two runs pinned to the same time steps and
# converged well past the reporting precision.  Left to itself the adaptive stepper takes
# different substeps in the two runs -- the DOF ordering alone is enough to change where
# each Newton iteration stops inside the tolerance band -- and the difference grows to a
# few parts in a thousand for reasons that have nothing to do with the aquifer.

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
  exit 1
fi

OPTIND=1
while getopts "i:f:r:a:t:c:e:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
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
${COMPARE_ECL_COMMAND} -t SMRY -a ${RESULT_PATH}/grid/${FILENAME} \
                                  ${RESULT_PATH}/aux/${FILENAME} \
                                  ${ABS_TOL} ${REL_TOL}
exit $?
