# Assembly transfer optimization: RSC_96 injection

The same case and solver settings from `compile_and_test.sh` were built and run
through `dockerrunscript.sh` on the ROCm image and AMD Radeon Pro W7900, with
one MPI rank and 16 CPU threads. The baseline simulator revision was
`4a116254e`. The source-cell change was subsequently checkpointed by the user
in `f53bb7f3d`; the remaining changes add resident true-IMPES CPR weights.

## Measured result

Three full baseline runs and three full optimized runs gave these timings
(seconds). Validation and profiling were disabled for timing runs.

| Run | Baseline assembly | Optimized assembly | Baseline simulation | Optimized simulation | Linear iterations, baseline / optimized |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 34.38 | 16.11 | 134.05 | 99.95 | 308 / 314 |
| 2 | 34.32 | 16.38 | 134.36 | 105.43 | 302 / 316 |
| 3 | 33.80 | 16.01 | 134.20 | 104.76 | 317 / 315 |
| Median | 34.32 | 16.11 | 134.20 | 104.76 | |

Assembly is **2.13 times faster** (53.1% less time). Total simulation takes
**21.9% less time**. All six runs completed 57 timesteps, 255 linearizations and
198 Newton iterations, with no wasted iterations and 57 relaxed CNV acceptances.

Median linear-solve time decreased from 70.95 to 49.30 s. Pre/post-step time
increased from 20.87 to 31.40 s because remaining CPU consumers now materialize
properties later. The total simulation improvement above includes this cost;
the result does not rely solely on moving work outside the assembly timer.

## Transfers and synchronization

1. **Thermal source evaluation downloaded every cell's intensive quantities
   (IQs).** The deck supplies injection temperature rather than explicit heat
   rate, so enthalpy evaluation asks for CPU IQs. Previously, requesting one
   source cell materialized all 903,168 cells, about 657.5 MB per cache slot.
   The new source accessor copies only the requested cell, overlays the same
   CPU fields and marks only that cell valid. Phase calculations reuse that
   entry. Normal cache invalidation handles subsequent property updates.
2. **True-IMPES CPR weights were calculated on the CPU.** This later consumer
   also materialized the whole IQ cache and uploaded three weights per cell.
   A new kernel evaluates the shared TPFA storage calculation from resident
   IQs, forms the transposed storage-derivative block and uses the existing GPU
   block solver and the same normalization to produce resident weights. It
   preserves the configured true-IMPES algorithm. Unsupported configurations
   and singular-block cases retain the CPU calculation and its diagnostics.
3. **Ordering is explicit.** A GPU stream event orders the weight kernel after
   property writes. A four-byte status download per weight update checks for
   singular blocks; it replaces the large property/weight round trip. Source
   cell copies still wait before their CPU consumers use the data.

The source-only change reduced assembly to 15.89 s but left all 200 bulk IQ
copies in the later CPR path, yielding no overall improvement. That measurement
identified why the CPR change was necessary.

Every final timing run recorded the same transfer totals:

| Counter | Baseline | Optimized |
| --- | ---: | ---: |
| Full IQ-cache downloads | 200 | 59 |
| Source-cell IQ downloads | 0 (covered by full copies) | 22,968 |
| Total IQ download bytes | 131,501,260,800 | 38,809,592,640 |
| Source-cell subset of IQ bytes | 0 | 16,720,704 |
| Primary-variable downloads | 198 | 57 |
| Primary-variable download bytes | 8,583,708,672 | 2,471,067,648 |
| Successful resident CPR-weight updates | 0 | 198 |

The IQ counter includes bulk and source-cell transfers; source counters are
subsets, not additional traffic. IQ download traffic decreased by 70.5%.
All 198 weight callbacks used the resident result, removing their per-update
weight uploads (4,291,854,336 bytes, derived from 903,168 cells × 3 doubles ×
198 updates). The existing one-time initialization upload remains.

Remaining transfers include bulk IQ and primary-state materialization at CPU
boundaries, the full residual download required by CPU convergence checks,
small source-delta uploads, and the device-to-device residual flattening copy.
The no-well GPU solver path already retained the Jacobian on the device.
`waitForAssembly()` already used a GPU stream event, not a CPU wait; it was not
replaced with a device-wide synchronization.

## Validation and numerical limits

The HIP build passed. A complete run with `OPM_GPU_VALIDATE_TRUEIMPES=1`
compared **536,481,792 scalar weights** against the existing CPU computation
on the same states: all 198 solves and all 903,168 cells. The largest absolute
difference was **7.216449660063518e-16**, below the diagnostic's 1e-10 threshold.
This checks the shared storage derivatives, geometric scaling, transposition,
small-system solution and normalization over the complete injection schedule.
The diagnostic intentionally downloads IQs and weights, so it is unset for
performance measurements.

The existing strict restart/summary comparator **does not pass even between
unchanged baseline runs**. No exact-output-equivalence claim is made. All
saved restart records were additionally compared for structure, finiteness,
maximum absolute difference and relative L2 error. The optimized runs show
similar variability to baseline repeats. Relative L2 differences below are
against the first baseline and cover all saved reports/cells:

| Field | Baseline repeats | Optimized runs |
| --- | ---: | ---: |
| Pressure | 1.36e-6–2.82e-6 | 2.16e-6–3.09e-6 |
| Capillary pressure | 1.66e-5–1.91e-5 | 1.91e-5–1.92e-5 |
| Gas saturation | 2.52e-3–3.53e-3 | 2.50e-3–2.52e-3 |
| Water saturation | 1.24e-4–1.73e-4 | 1.23e-4–1.24e-4 |
| Temperature | 7.24e-5–1.12e-4 | 7.20e-5–7.22e-5 |

Localized differences are much larger than the global norms: unchanged
baseline repeats differ by up to 0.20 saturation and 4.24 K; optimized runs
show up to 0.20 saturation and 3.09 K. These comparisons characterize existing
run variability; they are not a tighter numerical accuracy guarantee. The
same-state CPU/GPU weight comparison provides the direct validation of the new
weight computation. Composition fields are absent from the case's saved restart
output. CUDA and other physical models were not exercised by these measurements.

## Reproduction and evidence

`compile_and_test.sh` and `dockerrunscript.sh` have been restored byte-for-byte
to their original contents. The ordinary command now uses the optimized binary:

```sh
bash dockerrunscript.sh
```

To repeat the weight diagnostic, export `OPM_GPU_VALIDATE_TRUEIMPES=1` **inside
the container workflow**, before the simulator command in `compile_and_test.sh`.
Unset it again for timings. The optional diagnostic reports the maximum error
and throws on a nonfinite result or any absolute discrepancy above 1e-10.

Evidence directory: `/tmp/assembly_optimization_20260923`.

- `final_timings.json`: all six timings, iteration totals, counters and medians.
- `variability_baseline_2.json`, `_3.json`, and
  `variability_candidate_1.json` through `_3.json`: saved-field comparisons.
- `compare_variability.py` and `summarize_timings.py`: analysis scripts.
- Original workflow scripts, baseline/candidate executables and SHA-256 hashes.
- `baseline/` and `candidate_sparse_only/`: preserved initial output.

Other complete outputs are under `build/opm-simulators/bin/`:
`assembly_20260923_baseline_2`, `_3`, `assembly_20260923_candidate_1` through
`_3`, and `assembly_20260923_weight_validation`.

Console logs are `/tmp/assembly_baseline_20260923.log`,
`/tmp/assembly_baseline_repeats_20260923.log`,
`/tmp/assembly_weights_validation_retry_20260923.log`, and
`/tmp/assembly_final_candidate_runs_20260923.log`.

Validated and timed executable (preserved as `flow_gpu_final_candidate`) SHA-256:
`a6ca42ab552d3029171b8438e96938381fd672e33edf5278ac0b9cbfddf03404`.

A subsequent relink rebuilt only `moduleVersion.cpp.o`, which contains build
metadata. Implementation sources match the tested snapshot and simulator
libraries were not rebuilt. `final_manifest.json` records both executable hashes
and the checks of unchanged solver settings, restored scripts and zero wasted
iterations. The preserved executable is the exact artifact used for validation
and timing.
