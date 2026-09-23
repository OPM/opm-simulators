# Assembly transfer investigation (2026-09-23)

The current `compile_and_test.sh`, launched through `dockerrunscript.sh`, runs
RSC_96_nodisp_nodiff_injection with 16 threads and GPU properties, assembly and
linear solving. The starting simulator revision is
`4a116254e` (the full revision and original workflow scripts are saved under
`/tmp/assembly_optimization_20260923`).

## Transfers identified in the current code

- The no-well GPU solver path already retains the Jacobian on the device.
  The fallback still downloads and uploads it to support CPU contributions.
- Property kernels hand off to assembly through a GPU stream event wait;
  `waitForAssembly()` does not wait on the CPU.
- Source assembly uploads source indices, residual contributions and Jacobian
  contributions. These are small arrays, but use synchronous copies.
- Assembly downloads the complete residual. CPU convergence calculations
  still use every cell's residual, so this transfer is currently required.
- `flattenedResidual()` performs a device-to-device copy for the solver.
- Thermal SOURCE enthalpy calculation calls the CPU IQ accessor. Previously,
  one source cell caused the full device IQ slot to be downloaded and all CPU
  entries to be populated. The injection deck supplies temperatures, so the
  heat-rate shortcut does not avoid this path. Even phases with zero mass
  rate currently evaluate enthalpy.

## Implemented change

Thermal source evaluation now requests `intensiveQuantitiesForSource()`.
For a valid GPU property cache, the model downloads only the requested cell,
using the property stream to order the read after property writes. It waits
for this small transfer before overlaying the CPU IQ fields, then marks only
that cell valid. Further phase calculations for the same cell reuse the
CPU entry. Existing cache invalidation handles property updates and timestep
changes. Bulk consumers retain their existing whole-slot materialization.

The shared enthalpy calculation, source time index and arithmetic are preserved.
CPU models without the sparse accessor use their existing IQ accessor.
Transfer diagnostics include `source_iq_downloads` and
`source_iq_download_bytes`; these are subsets of total IQ download counters.

## Validation status

The unmodified baseline was compiled before applying the change. Its executable
is preserved at `/tmp/assembly_optimization_20260923/flow_gpu_baseline` and its
console log is `/tmp/assembly_baseline_20260923.log`.

The baseline completed with 57 timesteps, 255 linearizations, 198 Newton
iterations, 308 linear iterations and zero wasted iterations. Assembly took
34.38 s; simulation took 134.05 s. IQ diagnostics reported 200 downloads,
totalling 131,501,260,800 bytes. Saved output is under
`/tmp/assembly_optimization_20260923/baseline`.

The candidate build/run log is `/tmp/assembly_candidate_20260923.log`.
Candidate validation is in progress; no speedup is established yet.
