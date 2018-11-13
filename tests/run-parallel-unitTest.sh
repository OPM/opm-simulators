#!/bin/bash 
# This executes a unit test in parallel.
NP=$1 
BDIR=$2 
shift 2 
mpirun -np $NP $BDIR/bin/$@
