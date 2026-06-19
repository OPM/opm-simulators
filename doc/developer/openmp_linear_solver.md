# OpenMP threading in the linear solver

This branch makes the DILU/CPR linear-solver path usable with shared-memory
(OpenMP) parallelism, so a simulation can use intra-node threads instead of (or
in addition to) MPI domain decomposition. The thread count is taken from
`OMP_NUM_THREADS`, as elsewhere in Flow.

All additions are opt-in and leave the existing serial and MPI paths unchanged.

## What is threaded

| Component | Change |
|-----------|--------|
| SpMV operators (`WellOperators.hpp`) | The interior-row matrix-vector products in `GhostLastMatrixAdapter` and `WellModelGhostLastMatrixAdapter` run `#pragma omp parallel for` over output rows. Each `y[i]` is written by one thread and the matrix is read-only, so the result is bit-identical to the serial apply. Active automatically when built with OpenMP. |
| Scalar product (`ThreadedScalarProduct.hpp`) | `ThreadedSeqScalarProduct` threads the dot/norm reduction, gated behind a 50k block-count threshold so it does not add fork/join overhead on small systems. Used automatically by `FlexibleSolver`. |
| Smoother (`DILU2.hpp`) | `MultithreadDILU2`, a multicolor DILU that scales with threads (see below). Opt-in via `--linear-solver=dilu2`. |
| CPR pressure stage (`AmgclPreconditioner.hpp`) | AMGCL smoothed-aggregation AMG on the OpenMP backend, scalar (1x1) only. Opt-in as the pressure-stage preconditioner `"amgcl"`. Requires building with AMGCL (see below). |

## Choosing a smoother: `dilu` vs `dilu2`

`--linear-solver=dilu` uses the existing `MultithreadDILU`, whose level-set
(wavefront) schedule produces one barrier per level (~120 per apply on a typical
grid). That does not scale and regresses past a few threads.

`--linear-solver=dilu2` uses `MultithreadDILU2`, which reorders the unknowns by a
graph coloring: rows of one color are independent, so each triangular solve needs
only `#colors` parallel sweeps (a handful for grid graphs). The apply scales to
~3.7x at 8 threads where wavefront DILU regresses.

Trade-off: DILU2 factorizes the color-permuted matrix, so it is a slightly weaker
preconditioner (typically +~30% iterations), but it actually scales, and the
coloring -- hence the iteration count -- is fixed regardless of thread count, so
results are reproducible across thread counts.

## Building with AMGCL (for the `amgcl` pressure stage)

AMGCL is header-only and optional. Point CMake at an AMGCL clone:

```
cmake -DAMGCL_ROOT=/path/to/amgcl ...
```

CMake then defines `HAVE_AMGCL` and the `"amgcl"` preconditioner is registered for
scalar (1x1) systems. Without `-DAMGCL_ROOT` the entire AMGCL path compiles out and
nothing else is affected.

> Note: the fast numeric-only re-setup used by `update()` relies on a `rebuild()`
> entry point added by a small patch to AMGCL's `amgcl/amg.hpp`. That patch lives in
> the AMGCL clone, not in this repository, and should be upstreamed or carried as a
> tracked patch alongside the AMGCL dependency.

AMGCL is registered both in the serial factory (single-process shared-memory AMG)
and in the MPI factory (per-rank AMGCL wrapped as a restricted-additive-Schwarz
block preconditioner).

## Performance notes (single node)

On a memory-bandwidth-bound problem the threaded path scales with `OMP_NUM_THREADS`
up to the memory-bandwidth ceiling. On asymmetric CPUs (e.g. Apple performance +
efficiency cores) per-kernel barriers wait on the slowest core, so scaling can
plateau below the nominal thread count; uniform many-core nodes scale further. For
the CPR pressure stage specifically, very small pressure systems may not benefit
from more than ~4 threads.

For the trade-offs vs MPI and the rationale behind these choices, see the
threaded-linear-solver status notes accompanying this work.
