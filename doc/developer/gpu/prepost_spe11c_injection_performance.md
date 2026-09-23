# SPE11C injection pre/post processing

Measured on 2026-09-23 with ROCm 7, Radeon Pro W7900 (gfx1100), 16 threads,
and the `RSC_96_nodisp_nodiff_injection` case selected by
`compile_and_test.sh`. The case has 903168 active cells and runs to day 18250.
Baseline source: `28f367650`.

## Where the time went

Temporary wall-clock instrumentation reproduced the reported pre/post cost:
26.45 s, versus 26.43 s in the saved uninstrumented baseline. Instrumentation
has been removed from the final source.

| Profiled operation | Accumulated time |
| --- | ---: |
| Begin timestep | 2.33 s |
| Base end timestep, including well rate conversion and host materialization | 8.75 s |
| Summary evaluation | 12.62 s |
| Cell-output preparation, inside summary/output processing | 3.82 s |
| Separate fluid-in-place traversal, inside summary/output processing | 8.17 s |

The last two rows include initial output and overlap the summary row; they
must not be added to it. The remaining pre/post time includes time-level
advancement and cache copying.

Every fluid-in-place cell lookup acquired the GPU host-cache mutex, causing
contention in the OpenMP traversal. It also reread the large IQ cache after
the cell-output traversal had already loaded those properties. Well setup
and completion performed full-grid formation-factor and rate-conversion
scans even though this deck has SOURCE terms and no wells.

Some time charged to well completion was necessary GPU-to-host
materialization. Removing the rate-conversion scan moves that transfer to
the output boundary; it does not eliminate the whole 8.75 s. In an intermediate
measurement, output preparation increased by 6.69 s when that transfer moved.

## Changes

- Calculate fluid-in-place values in the existing cell-output traversal,
  using the intensive quantities already loaded into the element context.
  This removes a second grid pass and its contended per-cell cache lookups.
- Skip well formation-factor and rate-conversion scans when the global
  schedule contains no wells at the current report step and no actions.
  Schedules with shut wells, ranks without local wells, and actions that
  could create wells retain the existing calculations and collectives.
- Reuse the host IQ staging buffer without filling it with the prototype
  before each full GPU download. Newly allocated elements are still initialized.

No output frequency, solver settings, tolerances, physics, or cache-validity
protocol is changed.

## Timing results

| Metric | Saved baseline | Instrumented baseline | Final run 1 | Final run 2 |
| --- | ---: | ---: | ---: | ---: |
| Pre/post step | 26.43 s | 26.45 s | 19.03 s | 19.71 s |
| Simulation | 99.69 s | 98.83 s | 91.14 s | 92.71 s |
| Timesteps | 57 | 57 | 57 | 57 |
| Newton iterations | 198 | 198 | 198 | 198 |

The pre/post reduction is 25.4–28.0% against the saved baseline. These runs
use the original default GPU solver (Hypre). All builds and simulations run
through `bash dockerrunscript.sh`.

## Validation and local evidence

The original and final binaries were also run with `--linear-solver=dilu`
on the complete RSC_96 injection case. `compare_spe11c.py` reports zero
numerical differences for all 11 saved restart states (PCGW, PRESSURE, SGAS,
SWAT, TEMP) and all 104 physical summary samples. Complete convergence
histories match, including timestep retries: both runs report 136 timesteps,
820 linearizations, 716 Newton iterations, and 29133 linear iterations.
Composition fields are absent from this deck's restart output. Execution
speed summary keywords are excluded from numerical comparisons.

Hypre runs already exhibit nonzero restart and summary differences between
unchanged builds/runs, so those timing runs alone cannot establish numerical
equivalence. The separate DILU comparison avoids relying on that variability.

Two additional comparisons passed:

- Full RSC_16_nodisp_nodiff to day 730000 with DILU: all 27 restart states,
  all 2006 physical summary samples, and the complete convergence history
  match exactly. Pre/post time falls from 6.20 s to 3.35 s (46% lower), and
  simulation time from 19.79 s to 15.20 s. These are one run per binary.
- A 27-cell thermal CO2 injector/producer case exercises the well-present
  path and the existing CPU Newton fallback. All three restart files and
  both summary files are byte-identical, and convergence/retry histories
  match. Both runs report 33 timesteps and 354 Newton iterations.

The user-maintained `compile_and_test.sh` was restored byte-for-byte after
validation. No temporary timing instrumentation remains in the source.

Saved baselines, profiling logs, the validation driver, and the source diff
are under `/tmp/prepost_20260923`. Final run outputs are under
`build/prepost_20260923/candidate96_1` and `candidate96_2` in the devtainer
workspace. The same directory contains `baseline96_dilu`,
`candidate96_dilu`, `baseline16_long`, `candidate16_long`, `wells_baseline`,
`wells_candidate`, the original binary, the small well deck, and the additional
validation driver. JSON comparison results are in `/tmp/prepost_20260923`.
CUDA and MPI execution have not been tested.

For example, from the devtainer workspace:

```sh
python3 opm-simulators/doc/developer/gpu/compare_spe11c.py \
  build/prepost_20260923/baseline96_dilu \
  build/prepost_20260923/candidate96_dilu
```

Final candidate executable SHA-256:
`7d6bef6c3ec67720a74b27f60c79f2195e4f8c0876c2108aa60dbf685c89bb1b`.
