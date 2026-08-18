# Mixed-precision linear solvers
This folder contains mixed-precision building blocks for Krylov subspace methods
and a highly optimized mixed-precision implementation of ILU0 and DILU preconditioned
bicgstab. Hopefully, this will inspire the exploration of mixed-precision algorithms
in OPM.

The initial implementations are specialized for 3x3 block-sparse matrices due to their
importance in reservoir simulation. The original implementation only works in serial.
The parallel implementation wraps the optimized mixed-precision matrix-vector multiplication
and ILU0/DILU preconditioners into building blocks for the ISTL-based implementation of
bicgstab. This sacrifices some performance for improved modularity. Extending the work to
block-sparse matrices of arbitrary block size is work in progress.

The mixed-precision solver is selected by the command-line options `--linear-solver=mixed-ilu0`
or `--linear-solver=mixed-dilu`. The command-line option `--matrix-add-well-contributions=true`
must also be set as the mixed-precision solver operates directly on block-sparse matrices, not
on linear operators as other OPM solvers do. For convenience, a wrapper similar to the one
below can be used.

``` bash
OMP_NUM_THREADS=1 mpirun -np 1 --map-by numa --bind-to core build/bin/flow \
    --matrix-add-well-contributions=true \
    --linear-solver=mixed-ilu0 \
    --linear-solver-reduction=1e-3 \
    --linear-solver-max-iter=1024 \
    $@
```

To invoke the original serial implementation add the `legacy-` prefix to the mixed-precision
linear solver options.
