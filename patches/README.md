# AMGCL patch for the OpenMP CPR pressure stage

`AmgclPreconditioner` (see `opm/simulators/linalg/AmgclPreconditioner.hpp`) uses
AMGCL's `update()` -> `rebuild()` for cheap setup reuse between non-linear
iterations. AMGCL's stock `rebuild()` works but recomputes each coarse operator
with a full symbolic+numeric SpGEMM (`coarse_operator`), which is expensive and
barely threads.

`amgcl-numeric-galerkin-rebuild.patch` replaces that with a **numeric-only
Galerkin re-setup**: on the first build it records the coarse sparsity and a
"restriction map" (for each coarse nonzero, the list of `R_val*P_val` weights and
the index into the fine matrix values); subsequent `rebuild()` calls recompute the
coarse values only, parallel over coarse nonzeros (no atomics). The transfer
operators and aggregation are reused. This makes the re-setup ~9x cheaper and lets
it scale in OpenMP, which is what makes the AMGCL CPR pressure stage competitive.

The patch only touches `amgcl/amg.hpp` and hooks into the existing upstream
`rebuild()` / `allow_rebuild` machinery (it does not add public API). It is valid
for plain-Galerkin coarsenings (smoothed_aggregation, ruge_stuben) with a scalar
value_type, which is the configuration used here.

## Applying

AMGCL is a header-only dependency pointed to by `-DAMGCL_ROOT=<amgcl clone>`. From
that clone:

```
git apply /path/to/opm-simulators/patches/amgcl-numeric-galerkin-rebuild.patch
# or:  patch -p1 < /path/to/.../amgcl-numeric-galerkin-rebuild.patch
```

The patch header records the upstream AMGCL commit it was generated against.

## Status / upstreaming

This is carried as a patch rather than a vendored fork so the AMGCL dependency
stays a clean upstream checkout. It should be proposed upstream to AMGCL; until
then, building the `amgcl` preconditioner requires applying it. Building without
AMGCL (`-DAMGCL_ROOT` unset) compiles the whole path out, so this patch is not
needed for the rest of the branch.
