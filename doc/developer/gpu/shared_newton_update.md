# Shared black-oil Newton update

`BlackOilNewtonUpdate<TypeTag>::update()` contains the per-cell primary-variable
update previously implemented in `BlackOilNewtonMethod::updatePrimaryVariables_`.
The nonlinear method still runs on the CPU and owns iteration control, switch
history, and the count of cells whose primary-variable meanings changed.

The helper receives the problem, fluid-system instance, candidate/current primary
variables, correction, Newton parameters, and previous switch flag. It applies
the existing pressure and saturation limits, extension-specific bounds, phase
switching with oscillation suppression, and optional saturation projection. It
returns the new switch flag. The caller must supply a candidate carrying the
current cell's primary-variable meanings and metadata, as in the existing CPU
update path. Passing the same current and candidate object remains supported.

Primary-variable adaptation accepts an explicit fluid-system instance so future
device callers can use the non-static fluid-system view. The existing overload
continues to use the static CPU fluid system. Capillary-pressure evaluation uses
the existing device-callable black-oil fluid state for its saturation inputs.
Only models with a composition-switch variable instantiate oil/gas switching.

This extraction introduces no device Newton kernel, resident-state ownership,
new history API, floating-point compiler flags, or pressure-scaling behavior.
Those belong to subsequent GPU changes.

## Validation

The extraction is based on `origin/master` at `b28d5ccf1b`. Reproduction drivers
and logs are in the workspace's `testing/pr4-shared-newton-20261020/` directory.
All builds and runs use `dockerrunscript.sh`.

- Built CPU `flow_blackoil` and `flow_gaswater_energy`, HIP `flow_gpu` (gfx1100),
  and `test_blackoilprimaryvariables`.
- Passed the existing primary-variable test executable, including new checks for
  combined saturation/pressure limiting and switch-history suppression.
- Compared baseline and candidate CPU SPE1 summaries at zero absolute/relative
  tolerance, both with default pressure scaling and `--pressure-scale=100000`.
  This fixture does not request restart output; its restart files were absent.
- Compared short RSC_16 thermal source-only restart and physical summary output
  at zero absolute/relative tolerance. All comparisons use compareECL's standard
  summary blacklist. The source case includes injection followed by shutoff.
- Completed the HIP RSC_16 smoke run with GPU assembly and gpuISTL/DILU:
  54 linearizations, 27 Newton iterations and 37 linear iterations, no wasted
  iterations. The initial invocation omitted the GPU solver's required
  `--matrix-add-well-contributions=true`; retrying with that option passed.
- Verified that the extracted numerical limiting body matches the base after
  argument substitution and equivalent sign-function expansion.

No CUDA build, multi-rank run, device Newton kernel execution, or full extension
model regression suite was run. The HIP smoke checks the existing GPU assembly
path with CPU Newton updates; it does not validate a future resident Newton path.
