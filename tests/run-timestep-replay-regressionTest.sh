#!/bin/bash

# This runs a simulator, extracts the accepted substep end times from the
# generated INFOSTEP file, reruns the same case using the hardcoded timestep
# controller, and compares the two outputs.
#
# Note: When testing cases with timestep chops and cuts, passing --cpr-reuse-setup=0
# ensures bitwise identical Krylov solver trajectories between baseline (with chops)
# and hardcoded replay (without chops), isolating physical state rollback verification
# from heuristic preconditioner setup reuse cadences.

set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

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
  echo -e "\tOptional options:"
  echo -e "\t\t -d <path>     Unused, accepted for compatibility with other drivers"
  echo -e "\t\t -u <name>     Unused, accepted for compatibility with other drivers"
  echo -e "\t\t -y <option>   Additional option passed to replay simulation only"
  exit 1
fi

OPTIND=1
declare -a TEST_ARGS_REPLAY=()
while getopts "i:r:b:f:a:t:c:d:e:u:y:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    d) : ;;
    e) EXE_NAME=${OPTARG} ;;
    u) : ;;
    y) TEST_ARGS_REPLAY+=("${OPTARG}") ;;
  esac
done
shift $(($OPTIND-1))
if [ "${1:-}" = "--" ]; then
  shift
fi
declare -a TEST_ARGS=("$@")

BASELINE_PATH=${RESULT_PATH}/baseline
REPLAY_PATH=${RESULT_PATH}/replay
TIMESTEP_FILE=${RESULT_PATH}/${FILENAME}.timesteps
BASELINE_LOG=${RESULT_PATH}/${FILENAME}.baseline.log
REPLAY_LOG=${RESULT_PATH}/${FILENAME}.replay.log
BASELINE_INFOSTEP=${BASELINE_PATH}/${FILENAME}.INFOSTEP

resolve_simulator_binary() {
    # The test framework passes the full path to the simulator binary via -e.
    if [ -x "${EXE_NAME}" ]; then
        printf '%s\n' "${EXE_NAME}"
        return 0
    fi

    # Fall back to a flow_blackoil binary sitting next to a requested "flow", or dynamic flow binary.
    local bindir base
    bindir=$(dirname -- "${EXE_NAME}")
    base=$(basename -- "${EXE_NAME}")
    if [ "${base}" = "flow" ] && [ -x "${bindir}/flow_blackoil" ]; then
        printf '%s\n' "${bindir}/flow_blackoil"
        return 0
    fi
    if [ "${base}" = "flow_comp" ] && [ -x "${bindir}/flow_comp3_2p" ]; then
        printf '%s\n' "${bindir}/flow_comp3_2p"
        return 0
    fi
    if [ -x "${bindir}/flow" ]; then
        printf '%s\n' "${bindir}/flow"
        return 0
    fi

    echo "ERROR: Simulator binary not found: ${EXE_NAME}" >&2
    return 1
}

run_simulation() {
    local output_path=$1
    local log_path=$2
    shift 2
    local simulator_binary

    simulator_binary=$(resolve_simulator_binary) || return 1

    mkdir -p "${output_path}"
    "${simulator_binary}" "$@" --output-dir="${output_path}" "${INPUT_DATA_PATH}/${FILENAME}" > "${log_path}" 2>&1
    local status=$?
    if [ ${status} -ne 0 ]; then
        cat "${log_path}"
        return ${status}
    fi

    return 0
}

extract_timesteps() {
    local infostep_path=$1
    local log_path=$2
    local output_file=$3

    python3 "${SCRIPT_DIR}/extract_infostep_timesteps.py" "${infostep_path}" "${log_path}" "${output_file}"
}

rm -rf "${RESULT_PATH}"
mkdir -p "${RESULT_PATH}"

# Generate precise accepted substep end times from INFOSTEP.
run_simulation "${BASELINE_PATH}" "${BASELINE_LOG}" ${TEST_ARGS[@]+"${TEST_ARGS[@]}"} --output-extra-convergence-info=steps
test $? -eq 0 || exit 1

extract_timesteps "${BASELINE_INFOSTEP}" "${BASELINE_LOG}" "${TIMESTEP_FILE}"
test $? -eq 0 || exit 1

run_simulation "${REPLAY_PATH}" "${REPLAY_LOG}" ${TEST_ARGS[@]+"${TEST_ARGS[@]}"} ${TEST_ARGS_REPLAY[@]+"${TEST_ARGS_REPLAY[@]}"} --output-extra-convergence-info=steps --time-step-control=hardcoded --time-step-control-file-name="${TIMESTEP_FILE}"
test $? -eq 0 || exit 1

ecode=0
echo "=== Executing comparison for EGRID, INIT, UNRST, UNSMRY and RFT files ==="
"${COMPARE_ECL_COMMAND}" "${BASELINE_PATH}/${FILENAME}" "${REPLAY_PATH}/${FILENAME}" "${ABS_TOL}" "${REL_TOL}"
if [ $? -ne 0 ]
then
  ecode=1
  "${COMPARE_ECL_COMMAND}" -a "${BASELINE_PATH}/${FILENAME}" "${REPLAY_PATH}/${FILENAME}" "${ABS_TOL}" "${REL_TOL}"
fi

exit ${ecode}
