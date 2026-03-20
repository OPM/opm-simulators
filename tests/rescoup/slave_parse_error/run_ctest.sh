#!/bin/bash
# CTest wrapper for the reservoir coupling slave parse error test.
# Exits 0 on expected (currently unfixed) behavior, 1 on unexpected behavior.
#
# The slave deck has a deliberate parse error. The master should never
# complete successfully (exit 0). Expected outcomes depending on MPI:
#   - exit 1:   MPI runtime detected slave exit (OpenMPI 4.1.6)
#   - exit 9:   MPI runtime killed master when slave exited (custom MPICH v4.3.2)
#   - exit 124: master hung on MPI_Recv until timeout (OpenMPI 5.0.8, system MPICH v4.2.1)
#
# Any non-zero exit code is treated as PASS. Only exit 0 (master completed
# despite slave parse error) is treated as FAIL.
#
# Usage: run_ctest.sh <flow_binary> <mpi_launcher>

set -u

FLOW_BINARY="${1:?Usage: $0 <flow_binary> <mpi_launcher>}"
MPI_LAUNCHER="${2:?Usage: $0 <flow_binary> <mpi_launcher>}"

echo "Flow binary:  $FLOW_BINARY"
echo "MPI launcher: $MPI_LAUNCHER"
echo "Working dir:  $(pwd)"
echo ""

# Run with a timeout as a safety net. Some MPI implementations detect the
# slave exit immediately (exit 1 or 9), others hang until the timeout (exit 124).
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
else
    echo "PASS: Master did not complete (exit code $EXIT_CODE)"
    exit 0
fi
