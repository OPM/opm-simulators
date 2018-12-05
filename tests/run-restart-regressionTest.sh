#!/bin/bash

# This runs a simulator from start to end, then a restarted run of the
# simulator, before comparing the output from the two runs.  This is
# meant to track regressions in the restart support.

INPUT_DATA_PATH="${1}";       shift
RESULT_PATH="${1}";           shift
BINPATH="${1}";               shift
FILENAME="${1}";              shift
ABS_TOL="${1}";               shift
REL_TOL="${1}";               shift
COMPARE_ECL_COMMAND="${1}";   shift
OPM_PACK_COMMAND="${1}";      shift
PARALLEL="${1}";              shift
EXE_NAME="${1}";              shift
CASENAME="${1}";              shift
OPM_RSTFILE=$(echo "${1:-true}" | tr "[:upper:]" "[:lower:]");

if [ $# -gt 0 ]
then
    # Caller provided OPM_RSTFILE.  Shift away that setting to enable
    # passing all other parameters on to ${EXE_NAME}.
    shift
fi

BASE_NAME=`basename "${CASENAME}_RESTART.DATA"`

rm -Rf "${RESULT_PATH}"
mkdir -p "${RESULT_PATH}"

cd "${RESULT_PATH}"
if test "$PARALLEL" -eq 1
then
    CMD_PREFIX="mpirun -np 4 "
else
    CMD_PREFIX=""
fi

${CMD_PREFIX} "${BINPATH}/${EXE_NAME}" "${CASENAME}.DATA" \
    --enable-adaptive-time-stepping=false \
    --enable-opm-rst-file="${OPM_RSTFILE}" \
    --output-dir="${RESULT_PATH}" "$@"
test $? -eq 0 || exit 1

${OPM_PACK_COMMAND} -o "${BASE_NAME}" "${CASENAME}_RESTART.DATA"

    --enable-adaptive-time-stepping=false \
    --enable-opm-rst-file="${OPM_RSTFILE}" \
    --output-dir="${RESULT_PATH}" "$@"
test $? -eq 0 || exit 1

ecode=0
echo "=== Executing comparison for summary file ==="
${COMPARE_ECL_COMMAND} -R -t SMRY \
    "${RESULT_PATH}/${FILENAME}" \
    "${RESULT_PATH}/${FILENAME}_RESTART" \
    ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
    ecode=1
    ${COMPARE_ECL_COMMAND} -a -R -t SMRY \
        "${RESULT_PATH}/${FILENAME}" \
        "${RESULT_PATH}/${FILENAME}_RESTART" \
        ${ABS_TOL} ${REL_TOL}
fi

echo "=== Executing comparison for restart file ==="
${COMPARE_ECL_COMMAND} -l -t UNRST \
    "${RESULT_PATH}/${FILENAME}" \
    "${RESULT_PATH}/${FILENAME}_RESTART" \
    ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]
then
    ecode=1
    ${COMPARE_ECL_COMMAND} -a -l -t UNRST \
        "${RESULT_PATH}/${FILENAME}" \
        "${RESULT_PATH}/${FILENAME}_RESTART" \
        ${ABS_TOL} ${REL_TOL}
fi

exit $ecode
