# SPE11 pre/post cache lookup optimization

Measured on 2026-09-22 using the existing ROCm 7 / gfx1100 Docker workflow,
16 threads, and the complete RSC_16_nodisp_nodiff case (day 730000).
Source baseline: opm-simulators `41460f994`.

## Cause and change

`FIBlackOilModel::cachedIntensiveQuantities()` called
`ensureHostIntensiveQuantities()` for each cell read. That method acquired
the host-cache mutex and checked every cell's validity, including when the
cache was already materialized. A grid traversal therefore incurred O(N²)
validity checks. `intensiveQuantities()` also performed the same check before
calling `cachedIntensiveQuantities()`.

Individual lookups now check only the requested entry under the existing
mutex. An invalid entry still triggers full-slot materialization, and explicit
whole-slot consumers still check the complete cache. The duplicate check in
`intensiveQuantities()` is removed. No new validity state or invalidation
protocol is introduced, and cache publication remains protected against
concurrent output readers.

## Results

| Metric | Saved baseline | Optimized run 1 | Optimized run 2 |
| --- | ---: | ---: | ---: |
| Pre/post step | 68.69 s | 4.40 s | 4.47 s |
| Simulation | 84.02 s | 18.19 s | 17.90 s |
| Assembly | 4.08 s | 2.58 s | 2.57 s |
| Properties/update | 1.01 s | 1.00 s | 1.02 s |
| Linear solve | 8.47 s | 8.87 s | 8.56 s |

The baseline is the pre-existing completed run saved before rebuilding;
both optimized measurements were run during this optimization. Both builds
and simulations used `bash dockerrunscript.sh` from the devtainer directory,
with its existing `compile_and_test.sh` settings unchanged. No OPM unit tests
were created or run. CUDA was not built or benchmarked.

`compare_spe11c.py` reports zero numerical differences for all 27 restart
states (PRESSURE, PCGW, TEMP, SWAT, SGAS), all 2006 physical summary samples,
and identical complete convergence histories against the saved baseline for
both runs. Composition fields are absent from this deck's saved restart
output. All runs used 2006 timesteps, 4134 linearizations, 2128 Newton
iterations, and 4149 linear iterations. Bridge transfer counters match,
including 2130 IQ downloads and zero correction downloads.

Local evidence is preserved under `/tmp/gpu_prepost_20260922`:
`baseline16`, `candidate16_1`, `candidate16_2`, and the two Docker logs
`candidate1.log` and `candidate2.log`. For example, from the source directory:

```sh
python3 doc/developer/gpu/compare_spe11c.py /tmp/gpu_prepost_20260922/baseline16 /tmp/gpu_prepost_20260922/candidate16_2
```

The final executable SHA-256 is
`11357145fc7f3b68cac36483a4cfe6712c70689ddc4cae1cff5330942ca54c67`.
