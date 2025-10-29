# Mixed-precision linear solvers
This folder contains mixed-precision building blocks for Krylov subspace methods
and a highly optimized mixed-precision implementation of ILU0 preconditioned bicgstab.
Hopefully, this will inspire the exploration of mixed-precision algorithms in OPM.

The initial implementation is specialized for 3x3 block-sparse matrices due to their
importance in reservoir simulation and restricted to serial runs. Extending the work
to parallel runs is a matter of extending OPM's parallel infrastructure and should be
relatively straight-forward.

The mixed-precision solver is selected by the command-line option `--linear-solver=mixed-ilu0`.
The command-line option `--matrix-add-well-contributions=true` must also be set as the
mixed-precision solver operates directly on block-sparse matrices, not on linear operators as
other OPM solvers do. For convenience, a wrapper similar to the one below can be used.

``` bash
OMP_NUM_THREADS=1 mpirun -np 1 --map-by numa --bind-to core build/bin/flow \
    --matrix-add-well-contributions=true \
    --linear-solver=mixed-ilu0 \
    --linear-solver-reduction=1e-3 \
    --linear-solver-max-iter=1024 \
    $@
```
