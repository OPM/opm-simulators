#!/bin/bash
# Copyright 2026 Equinor
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

# CTest wrapper for the multi-slave reservoir coupling parse error test.
#
# The master spawns two slaves: RES-1 (slave1) has a deliberate parse error and
# fails; RES-2 (slave2) is valid and reports OK, then waits for the master's
# go-ahead. The master must detect RES-1's failure and cleanly abort the healthy
# RES-2 (which is mid-handshake) without deadlocking on a collective
# MPI_Comm_disconnect, then exit with code 1.
#
# Expected outcome (with the master-hang fix): exit code 1 on all MPI
# implementations.
#   - exit 1:   master detected the slave failure and aborted cleanly  -> PASS
#   - exit 0:   master completed despite a slave parse error           -> FAIL
#   - exit 124: master hung (e.g. the old multi-slave disconnect bug)  -> FAIL
#   - other:    unexpected                                             -> FAIL
#
# Usage: run_ctest.sh <flow_binary> <mpi_launcher>

set -u

FLOW_BINARY="${1:?Usage: $0 <flow_binary> <mpi_launcher>}"
MPI_LAUNCHER="${2:?Usage: $0 <flow_binary> <mpi_launcher>}"

echo "Flow binary:  $FLOW_BINARY"
echo "MPI launcher: $MPI_LAUNCHER"
echo "Working dir:  $(pwd)"
echo ""

# Run with a timeout as a safety net. With the fix the master exits quickly with
# code 1; a timeout (124) would indicate the master hung (regression).
timeout 20 \
    "$MPI_LAUNCHER" -np 1 "$FLOW_BINARY" \
    RC_MASTER.DATA \
    --parsing-strictness=low \
    --output-dir=. \
    2>&1

EXIT_CODE=$?

if [ $EXIT_CODE -eq 1 ]; then
    echo "PASS: Master detected slave failure and aborted cleanly (exit code 1)"
    exit 0
elif [ $EXIT_CODE -eq 0 ]; then
    echo "FAIL: Master completed successfully — slave parse error was ignored"
    exit 1
elif [ $EXIT_CODE -eq 124 ]; then
    echo "FAIL: Master hung until timeout (possible multi-slave disconnect deadlock)"
    exit 1
else
    echo "FAIL: Master exited with unexpected code $EXIT_CODE"
    exit 1
fi
