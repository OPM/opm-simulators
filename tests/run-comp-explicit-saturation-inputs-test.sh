#!/bin/bash

# Check explicit saturation input through Flow using small static decks.
set -euo pipefail

while getopts "e:c:i:r:" OPT; do
  case "${OPT}" in
    e) SIMULATOR=${OPTARG} ;;
    c) CHECKER=${OPTARG} ;;
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    *) exit 1 ;;
  esac
done

if [[ -z ${SIMULATOR:-} || -z ${CHECKER:-} || -z ${INPUT_DATA_PATH:-} || -z ${RESULT_PATH:-} ]]; then
  echo "Usage: $0 -e simulator -c checker -i input-directory -r results" >&2
  exit 1
fi

run_case()
{
  local name=$1 normalize=$2
  local case_path="${RESULT_PATH}/${name}"
  mkdir -p "${case_path}"
  if ! "${SIMULATOR}" "${INPUT_DATA_PATH}/${name}.DATA" --output-dir="${case_path}" \
      --normalize-explicit-swat="${normalize}" > "${case_path}/flow.log" 2>&1; then
    cat "${case_path}/flow.log"
    return 1
  fi
  (cd "${case_path}" && "${CHECKER}" "${name}")
}

expect_failure()
{
  local name=$1 expected=$2
  local case_path="${RESULT_PATH}/${name}"
  mkdir -p "${case_path}"
  if "${SIMULATOR}" "${INPUT_DATA_PATH}/${name}.DATA" --output-dir="${case_path}" \
      > "${case_path}/flow.log" 2>&1; then
    echo "${name}: invalid initial state was accepted" >&2
    return 1
  fi
  if ! grep -Fq "${expected}" "${case_path}/flow.log"; then
    cat "${case_path}/flow.log"
    echo "${name}: expected input validation error" >&2
    return 1
  fi
}

run_case COMP_EXPLICIT_SOIL_RTEMP true
run_case COMP_EXPLICIT_SOIL_SGAS true
run_case COMP_EXPLICIT_DRY_SOIL false

# This input normalization must also run when flash SWAT normalization is off.
run_case COMP_EXPLICIT_SATURATIONS_SCALED false
warning='Normalized explicit saturations in 9 cells; largest absolute deviation of their input sum from 1: 0.02'
if ! grep -Fq "${warning}" "${RESULT_PATH}/COMP_EXPLICIT_SATURATIONS_SCALED/flow.log"; then
  echo 'Missing count or maximum deviation in the normalization warning' >&2
  exit 1
fi

expect_failure COMP_EXPLICIT_SATURATIONS_OVERFULL \
  'Invalid explicit saturations in cell (1, 1, 1): The supplied saturations sum to more than 1'
expect_failure COMP_EXPLICIT_INVALID_TEMPERATURE \
  'Invalid temperature in cell (1, 1, 1): expected a finite, positive value'
