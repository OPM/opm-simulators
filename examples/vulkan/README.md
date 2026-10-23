# Experimental Vulkan linear solver (RFC)

This standalone prototype runs Flow's black-oil TPFA model with a Vulkan FP64
linear solver. It is provided for architecture feedback, **not as a merge-ready
production backend or a claim of general GPU acceleration**. The normal OPM
build does not include this directory or acquire Vulkan dependencies.

Full-history Norne currently **fails the output comparison** despite both
backends completing. See the failure evidence below; this is a known blocker.

The motivation is an alternative GPU compute path on an AMD Radeon 8060S using
Mesa RADV. CUDA/HIP support and gpu-ISTL already exist in OPM; this proposal asks
whether Vulkan support is useful and where it should fit with that work.
Only this GPU and the installed OPM 2026.04 release have been tested so far.
The prototype is carried on current master, but compiling against the installed
release does **not** validate integration with master's headers and libraries.

## Architecture and scope

- Flow retains parsing, black-oil physics, wells, Jacobian assembly, Newton,
  adaptive time stepping and ECLIPSE-format output.
- A separate TypeTag selects an `AbstractISTLSolver` adapter. The CPU reference
  executable uses the same model and the native OPM solver.
- Vulkan runs FP64 BiCGStab, sparse products, vector operations, reductions,
  approximate DILU substitutions and the application of a pressure AMG cycle.
- CPU work includes true-IMPES weights from OPM, conversion and row scaling,
  DILU factorization, pressure hierarchy/value updates, coarse inverse setup
  and a final true-residual verification against the original matrix.
- CPR uses unsmoothed aggregation, a Galerkin hierarchy, weighted Jacobi and
  a dense coarse inverse. This is **not identical** to the CPU reference AMG.
- The stopping criterion uses original equation units. There is no CPU linear
  solver fallback. Failed convergence is reported to Flow.

Supported scope: three-equation blocks, one MPI process, Newton, and explicit
well contributions in the matrix. MPI distribution, NLDD, other block sizes,
other hardware, and general ECLIPSE compatibility are not validated. A unique
hardware device must match `--vulkan-device`; an empty match selects only if
there is exactly one hardware device. CPU Vulkan devices are rejected. The
runtime requires `shaderFloat64`, host-visible coherent allocations and enough
storage-buffer capacity; this memory strategy is not tuned for discrete GPUs.

The `Reslab` namespace and `RESLAB_*` environment variables identify the original
experimental implementation. They are intentionally retained in this RFC to
keep the numerical path close to the measured prototype. Naming, configuration
plumbing, library boundaries and integration should be settled with maintainers.

## Build and numerical tests

Requires CMake >= 3.23, a C++20 compiler, Vulkan development files,
`glslangValidator`, Python 3, and optionally installed OPM development packages.
Python is only used to embed compiled SPIR-V in C++ headers.

```sh
cmake -S examples/vulkan -B build-vulkan -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DOPM_VULKAN_GPU_TESTS=ON
cmake --build build-vulkan -j 2
OPM_VULKAN_TEST_DEVICE=8060S ctest --test-dir build-vulkan --output-on-failure
```

`BUILD_OPM_ADAPTER=OFF` builds just the runtime and numerical tests without OPM.
GPU tests default to OFF so that a build-only environment does not try to use
hardware. When enabled, the six tests run serially; lack of a matching GPU is
a failure, not a silent software fallback. They cover all four preconditioners,
cooperative CPR and fine-pressure-fine CPR, with known solutions, anisotropy,
weak connections, zero-edge activation, scaling, invalid input and iteration
limits. There is not yet a GPU-free numerical test target.

The adapter's `OPM_VULKAN_PACKAGE_ASSERTIONS=ON` matches the tested OPM packages
(`WITH_NDEBUG=OFF`). This is an ABI requirement: `EnsureFinalized` changes layout
with `NDEBUG`. Set it OFF when using release libraries built with `NDEBUG`.
The temporary version shim identifies the tested 2026.04 adapter; replacing it
with the normal moduleVersion target belongs to production-build integration.

`Containerfile` captures the tested release-package environment. Its OPM base
image is pinned by digest and OPM development package versions are fixed;
other distribution packages are not snapshot-pinned. Record the resulting
image ID and driver version when reporting runs.

## Paired deck runs

Use an input deck accepted by OPM, with both SUMMARY and restart output enabled.
The example uses the same tolerances, well treatment and 16 threads for both
backends. Run sequentially from the deck directory, substituting absolute paths.

```sh
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1
/path/to/build-vulkan/flow_vulkan_cpu MODEL.DATA \
  --output-dir=/path/to/cpu-output --threads-per-process=16 \
  --linear-solver=cpr_trueimpes --linear-solver-reduction=1e-8 \
  --linear-solver-max-iter=2000 --matrix-add-well-contributions=true \
  --ecl-output-double-precision=true
RESLAB_VULKAN_COOPERATIVE=1 RESLAB_VULKAN_AGGREGATION_STRENGTH=0.25 \
  /path/to/build-vulkan/flow_vulkan MODEL.DATA \
  --output-dir=/path/to/gpu-output --threads-per-process=16 \
  --linear-solver=vulkan-cpr --vulkan-device=8060S \
  --linear-solver-reduction=1e-8 --linear-solver-max-iter=2000 \
  --matrix-add-well-contributions=true --ecl-output-double-precision=true
python3 examples/vulkan/compare_results.py /path/to/cpu-output /path/to/gpu-output
```

The comparison script requires `numpy` and `resdata`. It compares every requested
SUMMARY vector at matching report dates and PRESSURE/SWAT/SGAS/RS/RV in matching
restart reports. Absolute tolerance is 0.02 **or** relative tolerance 1e-5;
neither is relaxed for GPU. Six solver performance counters are recorded
separately: MLINEARS, MSUMLINS, NLINEARS, NLINSMAX, NLINSMIN and TCPU. Their
units, shapes and finite values are still checked. Linear tolerances of 1e-4
and 1e-6 failed earlier Norne output comparisons and are not used here.

The public Norne input used in the original experiment required removal of a
fourth slash after the three TUNING records and printing controls. Both backends
used the same prepared deck with normal parsing strictness. Unsupported-option
warnings remain relevant; this does not certify equality with ECLIPSE.

## Evidence and remaining work

Before packaging this RFC, Norne (44,431 active cells, first 56 days) passed
613 checks in every measured run. One warmup plus three timed repetitions gave:

| Host threads | CPU CPR median | Vulkan CPR median |
|---|---:|---:|
| 1 | 7.481 s | 7.230 s |
| 16 | 5.280 s | 5.275 s |

The individual samples and input/binary hashes are in [measurements.json](measurements.json).
The 16-thread samples overlap: **a tie**, not a demonstrated multicore speedup.
These include startup, setup, simulation and output. They describe the original
prototype, not fresh measurements of the reformatted RFC binaries. SPE1 was
slower on GPU (single four-thread run: CPU 1.22 s, GPU 1.82 s; 536 checks).
Longer simulations do not eliminate repeated solver setup; larger grids may
amortize dispatch overhead but can worsen memory costs and convergence.

The subsequent full-history Norne run (3,312 days, 16 threads, original prototype)
completed in both backends, but **349 of 3,275 comparison checks failed** with
the same thresholds, including SUMMARY vectors and restart RS. The cause has
not been isolated. The observed times (CPU 227.06 s, GPU 209.39 s) are **not a
valid speedup**: numerical comparison failed, and the host also had development
build activity during the run. [Failure evidence](norne-full-validation.json).

The packaged RFC runtime passed all six CTest GPU tests on Radeon 8060S. The
standalone comparison script also reproduced all 613 successful checks from
the preserved short Norne outputs.

`prepare_refined_spe1.py SOURCE_SPE1CASE2.DATA TARGET_DIRECTORY` creates a
120,000-cell, ten-day scaling experiment. It refines the original 10 x 10 x 3
grid to 100 x 100 x 12 while preserving reservoir dimensions, layer properties
and PVT. It remaps the well completions; well indices change with refinement.
It preserves the original dataset license notice and records changes and hashes.
This is synthetic scaling evidence, not another industrial field validation.

With the packaged binaries, SPE1 completed with 536 checks passing. The refined
case completed in both backends but failed 18 of 72 checks (including pressure
and gas saturation). Repeating with linear tolerance 1e-10 in **both** backends
still failed 18 checks; the physical comparison thresholds were unchanged.
Adaptive paths also differed (8 CPU versus 7 GPU time steps at 1e-8), so the
cause needs investigation rather than attributing it solely to linear accuracy.
The next diagnostic should distinguish time-step/nonlinear-path differences
from errors in the linear operator or adapter. [Packaged validation evidence](packaged-validation.json).

Before proposing production integration:

1. Discuss whether Vulkan belongs in gpu-ISTL, another adapter, or an external
   project, and whether its maintenance and GPU CI costs are acceptable.
2. Build and test against matching development versions of OPM's modules.
   A root master configuration with the available 2026.04 dependencies failed
   (including missing `opm_add_executable`); the main suite is not validated.
3. Diagnose the full-history Norne comparison failure and validate larger cases;
   repeat timings only after numerical comparisons pass.
4. Integrate configuration, ABI/version handling and tests into OPM conventions;
   review allocation/dispatch limits and add GPU-free tests where practical.
5. Evaluate reuse of existing preconditioner infrastructure rather than adding
   another long-term AMG implementation without agreement.
