# Bulk GPU transfer performance on RSC_96 injection

Measured on 2026-09-23 through `bash dockerrunscript.sh`, with ROCm 7,
Radeon Pro W7900 (gfx1100), one MPI rank and 16 threads. The case is
`RSC_96_nodisp_nodiff_injection`, as selected by `compile_and_test.sh`.
The baseline was rebuilt from clean simulator revision `a20c77bdd`.

The user selected the full GPU path, including property evaluation. Both
baseline and candidate use `--experimental-compute-properties-on-gpu=true`.
The final compile/test script enables this mode. All recorded simulator
parameters match between the six timing runs, except output directories.
Solver selection, tolerances, physics, output frequency and the deck are unchanged.

## Implementation

- Download the assembled residual directly into its existing CPU vector.
  Previously each download allocated and initialized a temporary vector, then
  copied it into the destination. The removed staging vector was 21,676,032
  bytes for this case. Compile-time block-size assertions and runtime cell-count
  checks make the existing contiguous-layout requirement explicit. The transfer
  still completes before CPU convergence evaluation consumes the residual.
- Register the reusable host intensive-quantity staging vector as pinned memory.
  Queue its download on the property-writing stream and wait once before the CPU
  field overlay, replacing the separate event wait and pageable/default-stream
  download. Registration is released before storage resize, reset or destruction.
  This reuses the existing 657,506,304-byte vector; it does not allocate another
  full host cache. The pinned registration persists between downloads.

The source-cell accessor, cache invalidation, CPU field overlay and all numerical
calculations are unchanged. Required property and residual download bytes remain.
All transfer counters match across the six timing runs.

## Repeated timings

Times are seconds. Baseline and candidate runs alternate, without profiling or
shadow checks. Each completes 57 timesteps, 255 linearizations and 198 Newton
iterations, with no wasted iterations.

| Pair | Baseline assembly | Candidate assembly | Baseline simulation | Candidate simulation | Linear iterations, baseline / candidate |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 15.79 | 15.34 | 91.17 | 91.38 | 303 / 315 |
| 2 | 16.24 | 15.37 | 92.37 | 90.95 | 314 / 301 |
| 3 | 15.87 | 15.42 | 92.21 | 91.68 | 309 / 311 |
| Median | 15.87 | 15.37 | 92.21 | 91.38 | |

Assembly improves in every pair, with **3.2% lower median assembly time**.
Median simulation time is 0.9% lower, but that difference is smaller than the
baseline run-to-run range, so this last change alone does not establish a large
additional total-runtime gain. Median pre/post time is nearly unchanged
(19.35 to 19.27 s); these measurements do not isolate a benefit from pinning alone.
The default Hypre solver exhibits the previously documented run variability.

For context, the earlier same-case baseline at `4a116254e` had median simulation
134.20 s, assembly 34.32 s, linear solve 70.95 s, properties/update 6.64 s and
pre/post 20.87 s. The three final candidates have medians 91.38, 15.37, 49.00,
6.44 and 19.27 s respectively. This is **31.9% lower simulation time** across the
combined assembly, solver, cache and transfer changes. These are successive
benchmarks, not an additional randomized comparison. The recorded parameter maps
match, and the earlier measurements and correctness limits are documented in
[assembly_transfer_optimization.md](assembly_transfer_optimization.md) and
[prepost_spe11c_injection_performance.md](prepost_spe11c_injection_performance.md).
The small properties/pre-post differences against that early baseline should
not be presented as separately established speedups.

## Validation

The HIP build and `TestProductionGpuDispatcherContract` pass. The focused test
compares downloaded GPU properties with CPU values and verifies that the host
cache is actually overwritten. Three complete numerical comparisons pass:

- RSC_96 injection with DILU: all 11 saved restart states (PCGW, PRESSURE, SGAS,
  SWAT and TEMP), all 104 physical summary samples and the complete timestep/retry
  history match exactly. Both runs have 136 timesteps, 820 linearizations,
  716 Newton iterations and 29,133 linear iterations, including the same retries.
- Long RSC_16 run to day 730000 with DILU: all 27 saved states, 2,006 physical
  summary samples and the complete convergence history match exactly.
- Thermal CO2 injector/producer case: all three restart files and both summary
  files are byte-identical, and the complete convergence/retry history matches.
  This exercises the well-present fallback path.

References are the validated pre/post candidate outputs under
`build/prepost_20260923`. Their archived C++ patch was checked byte-for-byte
against the `28f367650..a20c77bdd` implementation diff. These extra cases are
correctness checks, not additional paired performance measurements. Composition
fields are absent from the selected restart output; execution-speed summary
keywords are excluded from the numeric comparison.

DILU is used only for deterministic regression checks; the timing runs and final
compile/test script retain the original default GPU solver, Hypre. Unchanged
Hypre runs already differ numerically, as documented in the earlier reports.
The full-GPU compile/test script was restored and checked after validation.

CUDA and multi-rank execution have not been tested in this session.

## Rejected experiment

A gather kernel batched source-cell IQ downloads, reducing their count from
22,968 to 198 without reducing bytes. Three alternating pairs instead increased
median assembly from 16.05 to 16.58 s. That implementation and its extra test
changes were removed. It is not part of the final candidate.

## Local evidence

From the devtainer workspace, `build/bulk_transfer_20260923` contains the paired
`baseline_1` through `_3` and `candidate_1` through `_3` outputs, the validation
drivers, `timings.json`, and `summarize.py`. The exact rebuilt baseline executable
is `build/source_batch_20260923/flow_gpu_baseline`, SHA-256:
`f585acaa7e68b1a8540e64d91641e22bcc6f0629452e9dc45100868c1900964c`.

Candidate executable SHA-256:
`b0ce6ba9907e4d46715012ed48e7716b90ee2d7011a30ca36a7719a6f80113d9`.

The main validation console log is `/tmp/bulk_transfer_20260923_validation.log`;
the extra-case log is `/tmp/bulk_transfer_20260923_extra_validation.log`.
`comparison96_dilu.json`, `comparison16_long.json`, `comparison_wells.json` and
`final_manifest.json` are in the bulk-transfer evidence directory. The exact
candidate executable is preserved there as `flow_gpu_candidate`.
The rejected batch experiment is preserved separately in `build/source_batch_20260923/batch_timings.json` and
`/tmp/source_batch_20260923_batch_complete.log`. Its obsolete `baseline_1` run
used an older binary and is excluded; the accepted batch comparison used the
freshly rebuilt `baseline_current_1` through `_3` outputs.
