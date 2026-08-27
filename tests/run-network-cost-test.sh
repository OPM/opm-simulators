#!/bin/bash

# Runs a deck and checks what the network solve cost, rather than what it
# computed. The regression tests already compare the answer; this one catches
# the failure that leaves the answer intact and the run three times more
# expensive -- a network solve that falls back every step, or one well being
# re-solved because the network and the well model disagree about it.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory:"
  echo -e "\t\t -i <path>      Directory to read the deck from"
  echo -e "\t\t -f <filename>  Deck file name, without .DATA"
  echo -e "\t\t -r <path>      Directory to write results to"
  echo -e "\t\t -e <filename>  Simulator binary"
  echo -e "\tBounds (each optional; checked only when given):"
  echo -e "\t\t -F <n>         Maximum network fallbacks to the relaxed update"
  echo -e "\t\t -N <n>         Maximum total Newton iterations"
  echo -e "\t\t -W <n>         Maximum total well solves"
  echo -e "\t\t -S <percent>   Maximum share of well solves on any single well"
  exit 1
fi

OPTIND=1
MAX_FALLBACKS=-1
MAX_NEWTON=-1
MAX_WELL_SOLVES=-1
MAX_WELL_SHARE=-1
while getopts "i:r:f:e:F:N:W:S:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    F) MAX_FALLBACKS=${OPTARG} ;;
    N) MAX_NEWTON=${OPTARG} ;;
    W) MAX_WELL_SOLVES=${OPTARG} ;;
    S) MAX_WELL_SHARE=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
"${EXE_NAME}" ${TEST_ARGS} --output-dir=${RESULT_PATH} "${INPUT_DATA_PATH}/${FILENAME}.DATA" \
    > ${RESULT_PATH}/run.log 2>&1
if test $? -ne 0
then
  echo "FAIL: the simulator did not finish; tail of ${RESULT_PATH}/run.log:"
  tail -20 ${RESULT_PATH}/run.log
  exit 1
fi

PRT="${RESULT_PATH}/${FILENAME}.PRT"
DBG="${RESULT_PATH}/${FILENAME}.DBG"
for f in "${PRT}" "${DBG}"
do
  if test ! -f "${f}"
  then
    echo "FAIL: expected output ${f} was not written"
    exit 1
  fi
done

# The network gives up on a step and reverts to the relaxed fixed-point update.
FALLBACKS=$(grep -c 'simultaneously is not possible' "${DBG}")
# One line per well solve.
WELL_SOLVES=$(grep -c 'inner iterations' "${DBG}")
NEWTON=$(grep 'Newton its=' "${PRT}" \
         | awk -F'its=' '{split($2,a,","); s+=a[1]} END{printf "%d", s}')
# The well taking the largest share of the well solves, and its share.
read TOP_WELL TOP_COUNT <<< $(grep 'inner iterations' "${DBG}" \
         | awk '{print $2}' | sort | uniq -c | sort -rn | head -1 | awk '{print $2, $1}')
if test "${WELL_SOLVES}" -gt 0
then
  TOP_SHARE=$(( 100 * TOP_COUNT / WELL_SOLVES ))
else
  TOP_SHARE=0
fi

echo "network cost for ${FILENAME} ${TEST_ARGS}:"
echo "  network fallbacks : ${FALLBACKS}"
echo "  Newton iterations : ${NEWTON}"
echo "  well solves       : ${WELL_SOLVES}"
echo "  busiest well      : ${TOP_WELL} with ${TOP_COUNT} (${TOP_SHARE}%)"

STATUS=0
check() {
  # name value bound
  if test "$3" -ge 0 && test "$2" -gt "$3"
  then
    echo "FAIL: $1 is $2, above the bound of $3"
    STATUS=1
  fi
}
check "network fallbacks" "${FALLBACKS}" "${MAX_FALLBACKS}"
check "Newton iterations" "${NEWTON}" "${MAX_NEWTON}"
check "well solves" "${WELL_SOLVES}" "${MAX_WELL_SOLVES}"
check "the busiest well's share of well solves" "${TOP_SHARE}" "${MAX_WELL_SHARE}"
exit ${STATUS}
