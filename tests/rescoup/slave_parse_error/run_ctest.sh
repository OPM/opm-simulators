#!/bin/bash
# CTest wrapper for the reservoir coupling slave parse error test.
# Exits 0 on expected (currently unfixed) behavior, 1 on unexpected behavior.
#
# The slave deck has a deliberate parse error. The master should never
# complete successfully (exit 0). Expected outcomes depending on MPI:
#   - exit 9:   MPI runtime killed master when slave exited (custom MPICH)
#   - exit 124: master hung on MPI_Recv until timeout (OpenMPI, system MPICH)
#
# Usage: run_ctest.sh <flow_binary> <mpi_launcher>

set -u

FLOW_BINARY="${1:?Usage: $0 <flow_binary> <mpi_launcher>}"
MPI_LAUNCHER="${2:?Usage: $0 <flow_binary> <mpi_launcher>}"

echo "Flow binary:  $FLOW_BINARY"
echo "MPI launcher: $MPI_LAUNCHER"
echo "Working dir:  $(pwd)"
echo ""

# Always run with a timeout. Custom MPICH (v4.3.2) exits immediately with code 9;
# OpenMPI (v5.0.8) and system MPICH (v4.2.1) hang until the timeout fires (exit 124).
timeout 10 \
    "$MPI_LAUNCHER" -np 1 "$FLOW_BINARY" \
    RC_MASTER.DATA \
    --parsing-strictness=low \
    --output-dir=. \
    2>&1

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "FAIL: Master completed successfully — unexpected for slave parse error"
    exit 1
elif [ $EXIT_CODE -eq 9 ]; then
    echo "PASS: MPI runtime killed master on slave exit (exit code 9)"
    exit 0
elif [ $EXIT_CODE -eq 124 ]; then
    echo "PASS: Master hung on MPI_Recv as expected (timed out after 10s)"
    exit 0
else
    echo "FAIL: Unexpected exit code $EXIT_CODE"
    exit 1
fi
