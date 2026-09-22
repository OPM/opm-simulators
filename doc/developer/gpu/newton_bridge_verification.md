# Resident Newton bridge verification

All simulator builds and executions use the existing ROCm Docker workflow:

```sh
cd /home/tobias/opm/docker_images/devtainer
bash dockerrunscript.sh
```

The active image is `multitalentloes/opm-assembly-preliminary:rocm7`, targeting
`gfx1100`. Run variations belong in `compile_and_test.sh`, inside this workflow.
The required simulator arguments are:

```sh
./flow_gpu \
    --threads-per-process=16 \
    --newton-min-iterations=1 \
    --matrix-add-well-contributions=true \
    --linear-solver-accelerator=gpu \
    --linear-solver=dilu \
    --enable-storage-cache=false \
    --experimental-compute-properties-on-gpu=true \
    --output-dir=UNIQUE_OUTPUT_DIRECTORY \
    /home/tobias/opm/RSC_cases/SPE11C_cases/RSC_16/deck/RSC_16_nodisp_nodiff
```

## Fresh unmodified baseline, 2026-09-22

The baseline was built and completed before any implementation source edits.
All three source repositories were initially clean:

| Repository | Revision |
| --- | --- |
| opm-simulators | `651649dd2d5c42268018feb407d972021ac26f76` |
| opm-common | `25913b86b253d77ffc214d438b89fe5f20a4abb3` |
| opm-grid | `685ce0c4626c92104852293dbe315acc9b6687a4` |

Output: `../build/opm-simulators/bin/gpu_newton_20260922_baseline`.
The full console log, binary, workflow scripts, output copies and SHA-256
manifest are preserved in `/tmp/gpu_newton_verification_20260922`.
The baseline executable SHA-256 is
`405fc2bd0d65098cf5acb8ef8f9a84bdcc6bbc48e659b43692e265908b3bd76d`.

| Metric | Fresh baseline |
| --- | ---: |
| End day | 730000 |
| Timesteps | 2006 |
| Linearizations | 4134 |
| Newton iterations | 2128 |
| Linear iterations | 4149 |
| Wasted iterations | 0 |
| Simulation time | 138.00 s |
| Properties/update time | 34.33 s |
| Assembly time | 9.29 s |
| Linear solve time | 8.87 s |
| Pre/post step time | 84.74 s |

This single initial run is correctness evidence, not a repeated performance
measurement. Final speedup must use three complete enabled and disabled runs.

## Saved output comparison

From the simulator source directory:

```sh
python3 doc/developer/gpu/compare_spe11c.py BASELINE_DIRECTORY CANDIDATE_DIRECTORY
```

The standalone comparator reads Eclipse binary files with Python's standard
library. It compares pressure in Pa, temperature in K, saturation and any
available composition fields with relative tolerance `1e-6`, absolute floors
`1e-2 Pa`, `1e-6 K`, and `1e-8` for dimensionless quantities. It also compares
summary samples and reports differences in timestep and iteration history.
The fresh baseline has 27 saved states and 2006 summary samples. Its restart
fields include `PRESSURE`, `PCGW`, `TEMP`, `SWAT`, and `SGAS`; the deck does not
include composition in the default Eclipse-compatible restart output. Run
both lifecycle comparisons with `--enable-opm-rst-file=true` to include the
extended `RSW` and `RVW` fields. The comparator reports missing composition
explicitly and cannot establish equivalence from absent fields.

The comparison excludes execution timing summary fields (`ELAPSED`, `TCPU`,
`TCPUDAY`, `TCPUTS`, `TELAPLIN`), whose values should change between runs. It
retains physical times, all requested physical summaries and iteration counts.

## Shared-logic regression

The full bridge-disabled run completed successfully with
`--experimental-gpu-newton-update=false`, output
`../build/opm-simulators/bin/gpu_newton_20260922_disabled`. It produced exactly
identical saved physical fields, physical summary samples and complete timestep
and iteration history to the fresh baseline (maximum numeric error zero).
Simulation time was 136.50 s and properties/update time 33.71 s; convergence
totals remained 2006 timesteps, 4134 linearizations, 2128 Newton iterations and
4149 linear iterations. This establishes preservation of the CPU shared update
logic for the full case. It does not measure the resident path.

The container helper `run_spe11c_bridge.sh` provides `shadow`, `lifecycle`,
`rollback`, and `timing` modes, or `acceptance` to run them sequentially with
distinct output directories. Invoke it from `build/opm-simulators/bin` inside
the existing `compile_and_test.sh` workflow:

```sh
bash ../../../opm-simulators/doc/developer/gpu/run_spe11c_bridge.sh acceptance gpu_newton_20260922
```

Per-update shadow comparison, forced-rejection lifecycle validation, transfer
counters, repeated timing, and a ROCm trace are separate required checks.
Passing output comparison alone does not establish those acceptance criteria.

## Resident validation and numerical fix

The resident Newton kernel now compiles in its own translation unit with
floating-point contraction disabled (`-ffp-contract=off` for HIP,
`--fmad=false` for CUDA). Property kernels retain their original compiler
settings. The generated HIP compile command was checked to confirm that the
flag applies only to `GpuBlackoilNewtonUpdate.hip`.

With this change, `gpu_newton_20260922_nofma_shadow` completed all 2128 resident
Newton updates, all per-update shadow comparisons, and all 16 synthetic
phase-switch/stabilization batches. All 27 saved states, including pressure,
temperature, saturation, RSW and RVW, and all 2006 physical summary samples
were bit-identical to the bridge-disabled reference. Complete convergence
history also matched: 2006 timesteps, 4134 linearizations, 2128 Newton
iterations and 4149 linear iterations. This is stronger than the requested
tolerances; those tolerances were not changed. A later local rounding
difference of `5.68e-14` remained within the per-update tolerance and did not
change any saved output or convergence history.

The shadow run took 140.93 s, with properties/update time 35.53 s. These
validation timings include deliberate downloads and are not performance
measurements. Comparison results are in
`/tmp/gpu_newton_verification_20260922/nofma_shadow_vs_disabled.json`.

The numerical investigation is retained here to explain the compiler setting:

The first resident run completed all per-update shadow checks and all 16
synthetic phase-switch/stabilization batches. Saved lifecycle output comparison
nevertheless failed, so the acceptance sequence was stopped before rollback
and timing measurements. The failed outputs are retained as
`gpu_newton_20260922_shadow` and `gpu_newton_20260922_lifecycle_disabled`;
`gpu_newton_20260922_lifecycle_enabled` is an interrupted partial run.

The first summary pressure mismatch occurs at sample 1057, day 383615, with
349.426 Pa error against a 21.929 Pa tolerance. Reports 0–11 have bit-identical
saved pressure, temperature, RSW and RVW. Report 12 first exceeds saved-state
tolerances. The resident run has 2127 Newton iterations, versus 2128 for the
disabled run. The first local discrepancy occurred at the preceding phase
switch at day 383250: Rsw differed by `2.13e-14` (CPU
`32.519636577114412`, GPU `32.519636577114433`). Although that discrepancy
passed the local tolerance, its subsequent trajectory exceeded the saved-state
tolerance. Disabling contraction only for Newton update physics resolved that
failure without changing the property evaluator or relaxing tolerances.

## Optimized lifecycle and transfer counts

The complete optimized run `gpu_newton_20260922_nofma_lifecycle_enabled`
also produced identical saved fields, physical summaries and convergence
history to `gpu_newton_20260922_nofma_lifecycle_disabled`. Shadow validation
was off. Its 2128 successful device Newton updates performed no correction
downloads or primary-variable uploads within the resident update/property
path; the runtime checks enforce these invariants around every update.
The final timing executable also imports the CPU correction history once when
entering the resident path, preserving SOR semantics when execution changes
from CPU to GPU. The table includes that explicitly counted activation transfer.

| Counted category | Disabled calls / bytes | Enabled calls / bytes |
| --- | ---: | ---: |
| Primary-variable uploads | 2132 / 471564288 | 4 / 884736 |
| Primary-variable mirror downloads | 2006 / 443695104 | 4134 / 914374656 |
| Correction downloads | 2128 / 235339776 | 0 / 0 |
| Correction-history uploads on activation | 0 / 0 | 1 / 110592 |
| IQ downloads | 2130 / 7145349120 | 2130 / 7145349120 |
| Static-data uploads | 83 / 3952288 | 83 / 3952288 |
| Bridge allocations | 121 / 17955648 | 121 / 17955648 |
| Compatibility solver correction allocations | 2128 | 0 |

The four remaining primary-variable uploads occur during initialization.
Bridge buffers retain their initial allocations throughout resident updates.
The targeted path eliminates 706019328 transferred bytes per full run.
Required CPU mirrors add 470679552 download bytes relative to the disabled
path. Including the additional 110592-byte correction-history upload, the net
reduction across the counted bulk transfer categories is 235229184 bytes.
IQ downloads remain about 7.15 GB per run and are unchanged.
Assembly still reconstructs GPU parameter buffers each linearization; those
allocations are outside the persistent bridge and remain a separate cost.
Median total assembly time was 8.87 s disabled and 8.78 s enabled. This includes
matrix work as well as parameter reconstruction; it is not an allocation-only
measurement.

## Forced-rejection check: documented trajectory difference

The first paired rejection runs completed at day 730000 after injecting a
single rejection on the second timestep attempt. Both used 2008 timesteps,
4134 linearizations, 2127 accepted Newton iterations and 4044 linear iterations.
Their states were identical through report 11, but the report 12 comparison
failed. The first summary pressure failure was 114.441 Pa against a 21.930 Pa
tolerance. The first local Rsw difference was again only `2.13e-14` during a
phase switch at day 383250. The acceptance script stopped automatically before
timing measurements. Detailed evidence is preserved in
`/tmp/gpu_newton_verification_20260922/rollback_differences.json`.
Across saved report states, maximum absolute differences were 186.157 Pa for
pressure, 0.002808 K for temperature, `8.80e-5` for saturation, 0.0126534 for
Rsw and `4.46e-10` for Rvw. Rvw remained within its requested absolute floor;
the other listed fields exceeded their original tolerances somewhere.

A bounded follow-up diagnostic evaluated identical pressure and temperature
on both processors. The device `pow(T, 1.5)` and one `log` result differed from
the CPU by one ulp, propagating into the switched Rsw value. Density, molar
volume, input pressure and temperature were identical. This localized the
remaining discrepancy to CPU/GPU transcendental arithmetic rather than an
observed rollback-state mismatch. The original strict retry-output gate did
not pass. After reviewing this evidence, the user explicitly chose to retain
the existing physics and document the retry differences instead of changing
the phase-switch rules or introducing new transcendental implementations.
The original comparison tolerances remain unchanged.

## Reproducing final measurements

Inside the Docker workflow, from `build/opm-simulators/bin`:

```sh
bash ../../../opm-simulators/doc/developer/gpu/run_spe11c_bridge.sh timing-profile gpu_newton_20260922_final
```

This runs three disabled/enabled pairs with tracing and shadow checks off,
compares their physical output, writes a median timing report, then runs two
complete cases with ROCm tracing confined to seconds 5–7. It refuses to
overwrite existing output directories; choose a new prefix when repeating it.
The existing evidence uses the prefix shown above. The final tested executable
SHA-256 is `085ea11d310ca55c66e5edc0e6b59a1e76e3588d8e80a715797ba37d7da78649`.
Its tracked diff, new files and file hashes are preserved as
`integrated.final.*` under `/tmp/gpu_newton_verification_20260922`.

## Measured timings

All six final timing runs completed the full case with 2006 timesteps, 4134
linearizations, 2128 Newton iterations and 4149 linear iterations. Each
disabled/enabled pair had identical saved physical fields, physical summaries
and complete convergence history. Shadow checks and profiling were disabled.

| Repetition | Disabled simulation | Enabled simulation | Disabled properties/update | Enabled properties/update |
| --- | ---: | ---: | ---: | ---: |
| 1 | 136.68 s | 138.37 s | 33.83 s | 33.14 s |
| 2 | 138.05 s | 134.84 s | 34.45 s | 32.12 s |
| 3 | 135.98 s | 135.86 s | 33.67 s | 32.35 s |
| Median | 136.68 s | 135.86 s | 33.83 s | 32.35 s |

The measured median speedup is **1.006×** for simulation and **1.046×** for
properties/update. The 0.82-second simulation median difference is smaller
than the run-to-run spread, so these measurements do not establish a material
overall runtime improvement. The transfer elimination is deterministic and
verified by counters. Machine-readable timings and per-pair numerical
comparisons are preserved as `final_timings.json` and
`final_timing_comparison_1.json` through `_3.json` in the evidence directory.

## ROCm trace corroboration

Both profiled runs completed the full case and produced identical physical
output and convergence history. `rocprofv3` recorded HIP calls, memory copies
and kernel names only during seconds 5–7, after initialization. The disabled
window contains 28 property updates and 56 assembly kernels; the enabled
window contains 31 resident Newton/property updates and 62 assembly kernels.

| Correlated trace pattern | Disabled window | Enabled window |
| --- | ---: | ---: |
| Async host-to-device copy on property stream | 28 (one/update) | 0 |
| Async device-to-host copy on default solver stream | 28 (one/update) | 0 |
| Resident Newton kernels | 0 | 31 |
| All `hipMalloc` calls | 1344 (48/update) | 1457 (47/update) |
| All `hipFree` calls | 2240 (80/update) | 2449 (79/update) |

Each of the 28 disabled upload calls and 28 correction download calls is
followed by the property kernel. API correlation IDs, stream IDs and kernel
ordering distinguish them from the synchronous assembly/convergence copies
on the default stream. The ROCm CSV schema in this image omits copy sizes;
221184 bytes per primary-variable upload and 110592 bytes per correction
download come from the bridge counters. Thus the trace corroborates call
elimination without attributing unrelated transfers merely by size. Remaining
allocation calls belong outside the persistent bridge and are not claimed to
have been eliminated.

Trace CSV directories are
`../build/opm-simulators/bin/gpu_newton_20260922_final_profile_disabled_trace`
and `..._enabled_trace`. Reproduce their analysis with:

```sh
python3 doc/developer/gpu/summarize_spe11c_trace.py DISABLED_TRACE_DIRECTORY ENABLED_TRACE_DIRECTORY
```

The result is preserved as `final_trace_summary.json` in the evidence directory.
Verification used the ROCm/gfx1100 Docker configuration. CUDA source and build
registration are present, but CUDA was not built or run in this session. No
unstable OPM test suite was used as an acceptance gate.
