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

# CTest wrapper for the reservoir coupling test of unsupported target
# vectors (CTest target rc_slave_udq_target).
#
# The slave deck is valid, but a UDQ uses GOPRT for a slave group, which the
# slave cannot report: the master run sets that group's targets.  The slave
# stops while reading its deck, notifies the master via MPI, and the master
# exits with code 1.  The test also checks that the slave logged the reason,
# so that a failure for some other cause does not pass.
#
# For the known OpenMPI 5.x BTL TCP IPv6 issue with MPI_Comm_spawn, see
# ../slave_parse_error/run_ctest.sh.
#
# Usage: run_ctest.sh <flow_binary> <mpi_launcher>

set -u

FLOW_BINARY="${1:?Usage: $0 <flow_binary> <mpi_launcher>}"
MPI_LAUNCHER="${2:?Usage: $0 <flow_binary> <mpi_launcher>}"

echo "Flow binary:  $FLOW_BINARY"
echo "MPI launcher: $MPI_LAUNCHER"
echo "Working dir:  $(pwd)"
echo ""

rm -f RES-1*.log RC_SLAVE.PRT RC_SLAVE.DBG

# Run with a timeout as a safety net in case of regression (master hanging).
timeout 10 \
    "$MPI_LAUNCHER" -np 1 "$FLOW_BINARY" \
    RC_MASTER.DATA \
    --parsing-strictness=low \
    --output-dir=. \
    2>&1

EXIT_CODE=$?

if [ $EXIT_CODE -ne 1 ]; then
    echo "FAIL: Unexpected exit code $EXIT_CODE (expected 1)"
    exit 1
fi

EXPECTED="UDQ FUOPT (DEFINE in"
if ! grep -a -q -F "$EXPECTED" RES-1*.log RC_SLAVE.PRT 2>/dev/null; then
    echo "FAIL: Master exited with code 1, but the slave did not log"
    echo "      the unsupported target vector use ('$EXPECTED')"
    exit 1
fi

echo "PASS: Slave stopped on the unsupported target vector, master exited with code 1"
exit 0
