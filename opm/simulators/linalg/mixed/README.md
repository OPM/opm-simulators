# Mixed-precision linear solvers
This folder contains mixed-precision building blocks for Krylov subspace methods
for block-sparse linear systems of equations with highly optimized implementations
for select block-sizes. The mixed-precision  sparse matrix-vector multiplications
(SPMV) and preconditioners (ILU0/DILU) are combined to provide ILU0/DILU and
CPR+AMG preconditioned bicgstab algorithms. In the latter case, only the second
stage of the CPR algorithm is performed in mixed-precision. Moreover, the algorithms
leverage an improved scalar product implementation that takes advantage of the fact
that for parallel runs all ghost cells are sorted after local cells.

The current implementations work for both serial and parallel runs and for any block
size > 1. However, only block-sizes 2,3, and 4 benefit from hand-optimized implementations,
and suboptimal performance is expected for other block sizes. To run the simulator with
mixed-precision ILU0+BiCGSTAB, you can modify the wrapper script below to your liking
``` bash
OMP_NUM_THREADS=1 mpirun -np 1 --map-by numa --bind-to core build/bin/flow \
    --matrix-add-well-contributions=true \
    --linear-solver=mixed-ilu0 \
    --linear-solver-reduction=1e-3 \
    --linear-solver-max-iter=1024 \
    $@
```
Currently, a JSON specification file is required to activate mixed-precision CPR+AMG, i.e.
use the wrapper script below
```
OMP_NUM_THREADS=1 mpirun -np 1 --map-by numa --bind-to core build/bin/flow \
    --linear-solver=../mixed-cprw.json \
    $@
```
and modify the following `mixed-cprw.json` file to your liking
```
{
    "maxiter": "1024",
    "tol": "0.001",
    "verbosity": "0",
    "solver": "mixed-bicgstab",
    "preconditioner": {
        "type": "cprw",
        "use_well_weights": "false",
        "add_wells": "true",
        "weight_type": "trueimpes",
        "pre_smooth": "0",
        "post_smooth": "1",
        "finesmoother": {
            "type": "mixed-ilu0",
            "relaxation": "1"
        },
        "verbosity": "0",
        "coarsesolver": {
            "maxiter": "1",
            "tol": "0.10000000000000001",
            "solver": "loopsolver",
            "verbosity": "0",
            "preconditioner": {
                "type": "amg",
                "alpha": "0.33333333333300003",
                "relaxation": "1",
                "iterations": "1",
                "coarsenTarget": "1200",
                "pre_smooth": "1",
                "post_smooth": "1",
                "beta": "0",
                "smoother": "ilu0",
                "verbosity": "0",
                "maxlevel": "15",
                "skip_isolated": "0",
                "accumulate": "1",
                "prolongationdamping": "1",
                "maxdistance": "2",
                "maxconnectivity": "15",
                "maxaggsize": "6",
                "minaggsize": "4"
            }
        }
    }
}
```

The legacy mixed-precision implementation is still available. Unlike the current
implementation, it does not leverage the ISTL-framework and consequently only works in serial.
Moreover, only block sizes 3 and 4 are currently supported. To run the simulator with legacy
mixed-precision ILU0+BiCGSTAB, you can modify the wrapper script below to your liking
``` bash
OMP_NUM_THREADS=1 mpirun -np 1 --map-by numa --bind-to core build/bin/flow \
    --matrix-add-well-contributions=true \
    --linear-solver=legacy-mixed-ilu0 \
    --linear-solver-reduction=1e-3 \
    --linear-solver-max-iter=1024 \
    $@
```

Have fun!
