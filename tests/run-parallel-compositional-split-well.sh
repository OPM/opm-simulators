#!/bin/bash

# Copyright 2026 SINTEF Digital
#
# This file is part of the Open Porous Media project (OPM).
#
# OPM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OPM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OPM.  If not, see <http://www.gnu.org/licenses/>.

# Verify collective rejection of an active compositional well whose
# connections span MPI ranks.

set -u

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <mpi-launcher> <numproc-flag> <flow-binary> <source-deck> <result-dir>"
  exit 2
fi

MPI_LAUNCHER="$1"
MPI_NUMPROC_FLAG="$2"
FLOW_BINARY="$3"
SOURCE_DECK="$4"
RESULT_DIR="$5"

mkdir -p "${RESULT_DIR}"
SPLIT_DECK="${RESULT_DIR}/SIMPLE_COMP_SSHIFT_SPLIT_WELL.DATA"
RUN_LOG="${RESULT_DIR}/split-well.log"

# Add a second injector connection at the opposite end of the 30-cell row.
# The two-rank partition puts these connections on different owners.
if ! awk '
  { print }
  $1 == "INJ" && $2 == 1 && $3 == 1 && $4 == 1 && $5 == 1 &&
      $6 == "OPEN" && $7 == "2*" && $8 == "0.0151" && $9 == "/" {
    print "INJ 30 1 1 1 OPEN 2* 0.0151 /"
    ++insertions
  }
  END { if (insertions != 1) exit 1 }
' "${SOURCE_DECK}" > "${SPLIT_DECK}"; then
  echo "FAIL: Could not add exactly one split-well connection to ${SOURCE_DECK}"
  exit 1
fi

set +e
# Allow the partitioner to split the injector so the run reaches the
# compositional model's distributed-well check.
timeout 20 \
  "${MPI_LAUNCHER}" "${MPI_NUMPROC_FLAG}" 2 "${FLOW_BINARY}" \
  "${SPLIT_DECK}" --allow-distributed-wells=true \
  --output-dir="${RESULT_DIR}" \
  > "${RUN_LOG}" 2>&1
EXIT_CODE=$?
set -e

cat "${RUN_LOG}"

if [ "${EXIT_CODE}" -eq 0 ]; then
  echo "FAIL: A compositional well split across ranks was accepted"
  exit 1
fi

if [ "${EXIT_CODE}" -eq 124 ]; then
  echo "FAIL: Split-well rejection timed out"
  exit 1
fi

EXPECTED="Distributed compositional wells are not supported: well 'INJ' has connections on multiple MPI ranks"
if ! grep -Fq "${EXPECTED}" "${RUN_LOG}"; then
  echo "FAIL: Expected diagnostic was not emitted"
  exit 1
fi

echo "PASS: Split compositional well was rejected collectively"
