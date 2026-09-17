#!/bin/bash

# Runs an undersaturated deck and fails if free gas appears in the
# fluid-in-place report.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options>"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read deck from"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -f <filename> Deck file name"
  echo -e "\t\t -e <filename> Simulator binary to use"
  exit 1
fi

OPTIND=1
while getopts "i:r:f:e:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
  esac
done

rm -rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
"${EXE_NAME}" --output-dir=${RESULT_PATH} --enable-vtk-output=false \
  "${INPUT_DATA_PATH}/${FILENAME}.DATA" > /dev/null
test $? -eq 0 || exit 1

PRT=${RESULT_PATH}/${FILENAME}.PRT
test -f ${PRT} || { echo "No PRT file ${PRT}"; exit 1; }

# The fifth ':'-separated field holds free gas, dissolved gas and total.
awk -F: '/CURRENTLY IN PLACE/ { n++; split($5, gas, " "); if (gas[1] + 0 > 1) { bad = 1; print } }
         END { if (n == 0) { print "No fluid-in-place report"; exit 1 } exit bad }' ${PRT}
