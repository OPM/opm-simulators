#!/bin/bash
# CTest wrapper for the reservoir coupling slave parse error test.
#
# The slave deck has a deliberate parse error. The slave notifies the master
# of the failure via MPI, and the master exits with code 1 on all MPI
# implementations. Exit 0 (master completed despite slave parse error)
# is treated as FAIL.
#
# Known issue — OpenMPI 5.x BTL TCP IPv6 bug (observed with Ubuntu package):
#   On systems where a network interface has both IPv4 and IPv6 addresses
#   (common with WiFi), OpenMPI 5.x BTL TCP registers only the IPv6
#   address for the interface when built with --enable-ipv6. When
#   MPI_Comm_spawn is used, the spawned child may connect back to the
#   parent using an IPv4 source address that the parent doesn't recognize,
#   causing the connection to be rejected ("dropped inbound connection"
#   or "UNREACHABLE").
#
#   This bug is NOT present in OpenMPI 4.1.6. It has been observed with
#   the Ubuntu 25.10 system package (5.0.8-8ubuntu1) but could not be
#   reproduced with custom source builds of the same version, suggesting
#   additional runtime factors (e.g., --with-verbs, --with-libfabric)
#   affect connection routing.
#
#   Root cause: opal/mca/btl/tcp/btl_tcp_component.c mca_btl_tcp_create()
#   iterates opal_if_list and takes the FIRST address for each interface,
#   which is IPv6 when both families are present. The source code has a
#   comment acknowledging this limitation:
#   https://github.com/open-mpi/ompi/blob/v5.0.8/opal/mca/btl/tcp/btl_tcp_component.c#L501-L515
#
#   Workaround: disable IPv6 in BTL TCP address selection:
#     export OMPI_MCA_btl_tcp_disable_family=6
#   or persistently in ~/.openmpi/mca-params.conf:
#     btl_tcp_disable_family = 6
#   Note: the parameter uses literal 4/6, NOT AF_INET (2) / AF_INET6 (10).
#
# Usage: run_ctest.sh <flow_binary> <mpi_launcher>

set -u

FLOW_BINARY="${1:?Usage: $0 <flow_binary> <mpi_launcher>}"
MPI_LAUNCHER="${2:?Usage: $0 <flow_binary> <mpi_launcher>}"

echo "Flow binary:  $FLOW_BINARY"
echo "MPI launcher: $MPI_LAUNCHER"
echo "Working dir:  $(pwd)"
echo ""

# Run with a timeout as a safety net in case of regression (master hanging).
timeout 10 \
    "$MPI_LAUNCHER" -np 1 "$FLOW_BINARY" \
    RC_MASTER.DATA \
    --parsing-strictness=low \
    --output-dir=. \
    2>&1

EXIT_CODE=$?

if [ $EXIT_CODE -eq 1 ]; then
    echo "PASS: Master detected slave failure (exit code 1)"
    exit 0
elif [ $EXIT_CODE -eq 0 ]; then
    echo "FAIL: Master completed successfully — unexpected for slave parse error"
    exit 1
else
    echo "FAIL: Unexpected exit code $EXIT_CODE (expected 1)"
    exit 1
fi
