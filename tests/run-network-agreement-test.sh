#!/bin/bash

# Runs one deck twice -- once with the legacy relaxed network update, once with
# the simultaneous solve -- and compares the summary vectors.
#
# The simultaneous solve is an alternative way to reach the same balance, not a
# change of physics, so on a deck where the relaxed update converges the two
# must agree. That is checkable without any external reference, which matters
# because most of these decks have none.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [flags for the simultaneous run only]"
  echo -e "\t\t -i <path>      Directory to read the deck from"
  echo -e "\t\t -f <filename>  Deck file name, without .DATA"
  echo -e "\t\t -r <path>      Directory to write results to"
  echo -e "\t\t -e <filename>  Simulator binary"
  echo -e "\t\t -c <filename>  compareECL binary"
  echo -e "\t\t -a <tol>       Absolute tolerance"
  echo -e "\t\t -t <tol>       Relative tolerance"
  echo -e "\t\t -b <flags>     Flags given to BOTH runs"
  exit 1
fi

OPTIND=1
BOTH_ARGS=""
while getopts "i:r:f:e:c:a:t:b:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    c) COMPARE_ECL=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    b) BOTH_ARGS=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
NEWTON_ARGS="$@"

run() {  # subdir, extra flags
  local dir="${RESULT_PATH}/$1"; shift
  mkdir -p "${dir}"
  "${EXE_NAME}" ${BOTH_ARGS} "$@" --output-dir="${dir}" \
      "${INPUT_DATA_PATH}/${FILENAME}.DATA" > "${dir}/run.log" 2>&1
  if test $? -ne 0
  then
    echo "FAIL: the $1 run did not finish; tail of ${dir}/run.log:"
    tail -20 "${dir}/run.log"
    exit 1
  fi
}

run relaxed
run simultaneous ${NEWTON_ARGS}

DBG="${RESULT_PATH}/simultaneous/${FILENAME}.DBG"
SOLVED=$(grep -c 'solved the production network' "${DBG}")
REFUSED=$(grep -c 'simultaneously is not possible' "${DBG}")
echo "${FILENAME}: the simultaneous solve ran ${SOLVED} times, handed back ${REFUSED} times"
if test "${SOLVED}" -eq 0
then
  # Without this the test would pass by comparing the relaxed update with
  # itself, which is what happens when the solve declines every tree.
  echo "FAIL: the simultaneous solve never ran, so this compares nothing"
  exit 1
fi

# -d: compare at report steps only. The two paths place their ministeps
# differently -- on one deck that alone made the summary vectors different
# lengths -- and where a ministep falls is not part of the answer.
"${COMPARE_ECL}" -t SMRY -d -a -y \
    "${RESULT_PATH}/relaxed/${FILENAME}" \
    "${RESULT_PATH}/simultaneous/${FILENAME}" \
    "${ABS_TOL}" "${REL_TOL}"
exit $?
