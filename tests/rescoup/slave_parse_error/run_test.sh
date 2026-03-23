#!/bin/bash
# Manual test runner: Reservoir coupling slave deck parse error
#
# For automated testing use run_ctest.sh (registered as CTest target
# rc_slave_parsing_err). This script is kept for manual inspection: it
# prints human-readable output including the slave log file, making it
# easier to diagnose behavior interactively.
#
# This test runs a reservoir coupling master that spawns a slave process
# whose DATA file has a deliberate parse error (missing INCLUDE file).
#
# For expected behavior per MPI implementation and known OpenMPI 5.x issues,
# see the comments in run_ctest.sh.
#
# Usage:
#   run_test.sh <flow_binary_path> [mpi_launcher] [timeout_seconds]
#
# Examples:
#   ./run_test.sh /path/to/build/bin/flow
#   ./run_test.sh /path/to/build/bin/flow mpiexec
#   ./run_test.sh /path/to/build/bin/flow mpirun 60

set -u

FLOW_BINARY="${1:?Usage: $0 <flow_binary_path> [mpi_launcher] [timeout_seconds]}"
MPI_LAUNCHER="${2:-mpirun}"
TIMEOUT_SECONDS="${3:-30}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Resolve to absolute path before cd'ing to the test data directory
FLOW_BINARY="$(realpath "$FLOW_BINARY" 2>/dev/null || readlink -f "$FLOW_BINARY")"

if [ ! -x "$FLOW_BINARY" ]; then
    echo "ERROR: Flow binary not found or not executable: $FLOW_BINARY"
    exit 2
fi

if ! command -v "$MPI_LAUNCHER" >/dev/null 2>&1; then
    echo "ERROR: MPI launcher not found: $MPI_LAUNCHER"
    exit 2
fi

# Run from the master data file directory so relative paths in the deck work
cd "$SCRIPT_DIR" || exit 2

# Clean up any previous output files
rm -f RC_MASTER.INFOSTEP RC_MASTER.UNSMRY RC_MASTER.SMSPEC
rm -f RES-1*.log
rm -f *.PRT *.DBG

echo "=== Reservoir Coupling: Slave Parse Error Test ==="
echo "Flow binary:   $FLOW_BINARY"
echo "MPI launcher:  $MPI_LAUNCHER"
echo "Timeout:       ${TIMEOUT_SECONDS}s"
echo "Working dir:   $(pwd)"
echo ""
echo "--- Running test ---"

timeout "$TIMEOUT_SECONDS" \
    "$MPI_LAUNCHER" -np 1 "$FLOW_BINARY" \
    RC_MASTER.DATA \
    --parsing-strictness=low \
    --output-dir=. \
    2>&1

EXIT_CODE=$?

echo ""
echo "=== Results ==="
if [ $EXIT_CODE -eq 124 ]; then
    echo "RESULT: TIMEOUT (exit code 124)"
    echo "  Master hung waiting for slave — this is the expected current behavior."
    echo "  The slave failed to parse its deck and exited without sending"
    echo "  initial data to the master."
elif [ $EXIT_CODE -eq 0 ]; then
    echo "RESULT: SUCCESS (exit code 0) — unexpected, master should not succeed"
elif [ $EXIT_CODE -eq 137 ]; then
    echo "RESULT: KILLED (exit code 137, SIGKILL)"
    echo "  Process was killed, possibly by timeout or MPI runtime."
elif [ $EXIT_CODE -eq 1 ]; then
    echo "RESULT: ERROR (exit code 1)"
    echo "  Master exited with error — may indicate slave failure was detected."
else
    echo "RESULT: Exit code $EXIT_CODE"
fi

# Check slave log file if it exists
echo ""
echo "=== Slave Log Files ==="
FOUND_LOG=0
for logfile in RES-1*.log; do
    if [ -f "$logfile" ]; then
        FOUND_LOG=1
        echo "--- $logfile ---"
        cat "$logfile"
        echo ""
    fi
done
if [ $FOUND_LOG -eq 0 ]; then
    echo "(no slave log files found — slave may have exited before output redirect)"
fi

exit $EXIT_CODE
