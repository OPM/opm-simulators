# Hysteresis configuration, history, and calculated data in Flow

## Scope and source snapshot

This is a focused follow-up to [Flow simulation-time data flow](flow-simulation-data-flow.md).
It concerns black-oil saturation-function hysteresis during simulation, including
Carlson, Killough, and WAG behavior. It assesses a smaller hysteresis snapshot
as an initial refactoring target, then considers moving the authoritative
history out of material-law parameters. No implementation change or numerical
comparison was made for this analysis.

Inspected source revisions: `opm-simulators` `14daa38a6131` and `opm-common`
`7091cb2018bf`. The previous data-flow document describes an older revision;
in particular, optional timestep snapshot/restore is now wired into Flow.

## Assessment

The proposed three-way split fits the implementation, with one qualification:
“calculated” does not mean “can be reconstructed from the current saturation.”
Some WAG quantities are calculated from *earlier scanning curves*. Their
calculation is historical and its inputs may no longer exist. They must remain
in the authoritative state unless a smaller set of sufficient historical
values is identified and tested.

| Category | Examples in `EclHysteresisTwoPhaseLawParams` | How to treat it |
| --- | --- | --- |
| Configuration, fixed during an ordinary simulation step | `EclHysteresisConfig`, WAG configuration, drainage/imbibition curve parameters, endpoint scaling, SATNUM/IMBNUM selection, `Sncrd_`, `Sncri_`, `C_`, `Cw_`, `KrndMax_`, `Krwi_snmax_`, `pcmaxd_`, `pcmaxi_`, `Swco_` | Construct from deck/region/cell data; share or reference where possible. “Static” here means fixed over a step, not necessarily immutable for every use of the manager. |
| Evolving history that distinguishes two paths ending at the same current saturation | `pcSwMdc_`, `pcSwMic_`, `krnSwMdc_`, `krwSwMdc_`, `initialImb_`, WAG branch/reversal saturation values and cycle status | Preserve across timesteps, rollback, restart, and state transfer. |
| Values calculated from configuration and sufficient history | Carlson shift `deltaSwImbKrn_`; Killough `Sncrt_`, `Swcrt_`, `KrndHy_`, `KrwdHy_`, `Krwd_sncrt_`; selected WAG curve evaluations | Recompute when the underlying history changes, or calculate at use sites after measuring cost. Do not independently serialize/copy if equivalence is established. |

This classification is about *information*, rather than whether a member is
currently assigned in `update()` or `updateDynamicParams_()`. A cached value
can be expensive enough to retain for evaluation while still being absent
from the rollback snapshot. Conversely, a value written by a calculation can
be essential history if that calculation consumed an earlier curve that has
since been replaced.

## Present ownership and lifecycle

The per-cell `MaterialLawParams` held by `EclMaterialLawManager` contain two
phase hysteresis parameter objects. Those objects hold the curves and their
configuration together with mutable extrema, reversal points, branch status,
and derived coefficients. Three-phase approaches normally have gas-oil and
oil-water histories; a two-phase approach selects one of gas-oil, oil-water,
or gas-water. Directional relperm/imbibition can add X/Y/Z parameter sets.

```text
accepted/restored reservoir solution
    -> current cell fluid state
    -> FlowProblem::updateHysteresis_() at timestep entry
    -> manager updates base and possibly X/Y/Z material-law parameters
    -> Flow invalidates and recomputes current intensive quantities
    -> relperm and capillary-pressure evaluations read the changed parameters
```

`FlowProblem::updateHysteresis_()` deliberately visits all local elements,
including overlap/ghost elements, to keep parallel copies synchronized. Its
return value currently says hysteresis is enabled, not whether any cell
changed, so the current intensive quantities are invalidated even if no
history value actually changed. The manager's per-cell `updateHysteresis()`
does return a change flag, but the Flow wrapper discards it. This is a
separate opportunity from shrinking the snapshot.

Rollback now has two modes:

- With `EnableStateRollback=false` (the default), failed-step handling restores
  the reservoir primary variables; hysteresis history is not copied back
  explicitly. At retry entry, history is advanced again from the restored
  solution. This relies on the update behavior and call order.
- With `EnableStateRollback=true`, `FlowProblem::beginTimeStep()` captures a
  snapshot before updating explicit history, and `updateFailed()` restores it.
  The manager captures base and directional hysteresis state. The same option
  also captures other explicit histories in `FlowProblem` and affects Newton
  state rollback. The manager snapshot is a timestep-only object and is not
  restart-serialized.

The manager snapshot holds arrays of `EclHysteresisDynamicState` for each
applicable two-phase system and direction. Each element currently copies 23
scalar fields, one integer, and three booleans. It largely mirrors the mutable
members of the parameter object rather than a minimal history representation.
The number of arrays depends on phase approach and directional settings, so
even a modest reduction per element can matter for a large grid.

Two other persistence paths are distinct from this in-memory snapshot:

- The material-law manager's `serializeOp()` walks the base per-cell
  `MaterialLawParams`; each two-phase parameter `serializeOp()` lists nine
  dynamic members. For example, it includes `pcSwMic_` but not `pcSwMdc_`,
  and it omits the WAG branch state. The manager method does not itself walk
  directional parameter arrays. A snapshot redesign must audit this path
  separately; the in-memory snapshot is more complete than the current
  internal serialization list. This is a coverage difference that warrants
  a restart test, not by itself proof of a failed restart.
- ECL restart/output uses six familiar extrema for the oil-water and gas-oil
  systems and reconstructs parameter state through `set*HysteresisParams()`
  and `update()`. Those six values are an interchange format, not a proof that
  they suffice to reconstruct WAG scanning history or every internal cache.

## Field-by-field reduction candidate

The following table refers to the fields of
`EclHysteresisDynamicState<Scalar>` in
`EclHysteresisTwoPhaseLawParams.hpp`. “Drop” means drop from the **snapshot**
after adding a rehydration step; it does not require removing a cached member
from the material-law parameter object.

| Fields | Proposed status | Reason or condition |
| --- | --- | --- |
| `krnSwMdc`, `krwSwMdc` | Keep | Historical extrema; `krwSwMdc` currently also feeds output even where it does not alter a relperm branch. |
| `pcSwMdc`, `pcSwMic`, `initialImb` | Keep | Capillary-pressure reversal and initial imbibition history. `pcSwMic` alone does not determine `pcSwMdc`. |
| `deltaSwImbKrn` | Drop for Carlson after validation | Inverse imbibition curve applied to drainage relperm at `krnSwMdc`, minus that saturation. Only relevant to the Carlson branch. |
| `Sncrt`, `Swcrt` | Drop for ordinary Killough after validation | Explicit formulas in `updateDynamicParams_()` use historic extremum, fixed endpoints, and `C_`/`Cw_`. WAG's separate `SncrtWAG` needs different treatment. |
| `KrndHy`, `KrwdHy`, `Krwd_sncrt` | Drop after validation | Direct drainage-curve evaluations at `krnSwMdc` or at the reconstructed `Sncrt`. Respect model 4 and wetting-fix branches. |
| `isDrain`, `krnSwWAG`, `krnSwDrainRevert`, `krnSwDrainStart`, `krnSwDrainStartNxt`, `krnSwImbStart`, `swatImbStart`, `swatImbStartNxt` | Keep initially for WAG | Current branch and reversal coordinates are path-dependent. Some individual entries may later be eliminated by an invariant proof. |
| `nState` | Keep initially for WAG | It records cycle progression. Evaluation distinguishes cycle 1, cycle 2, and later cycles; the update increments the count. A saturated 1/2/3+ representation is plausible but needs a test of all uses. |
| `SncrtWAG`, `cTransf`, `krnImbStart`, `krnImbStartNxt` | Keep initially for WAG | These are calculated, but the formulas can depend on the previous scanning branch or only run at a particular transition. Recomputing from today's extrema is not generally equivalent. |
| `krnDrainStart`, `krnDrainStartNxt`, `krnSwImbStart` | Candidate second pass for WAG | They appear derivable from retained reversal coordinates and drainage-curve evaluations/inversion. Prove validity at *every* WAG branch transition and after restore before dropping. |
| `wasDrain` | Candidate second pass for WAG | `update()` assigns it from the prior `isDrain_` immediately before `updateDynamicParams_()` uses it. It may be reducible to a local transition variable, but capture/restore and equality behavior must be checked. |

The table intentionally overlaps some WAG entries: “keep initially” is the
safe first pass; “candidate second pass” records a possible later reduction.
There is no justified minimal WAG snapshot yet. A source-level dependency
argument is necessary but insufficient because sentinels, branch transitions,
and exact numerical equivalence also matter.

For non-WAG paths, a first useful snapshot would have the four historic
saturation values in the first two rows (`krnSwMdc`, `krwSwMdc`, `pcSwMdc`,
`pcSwMic`) and two status flags (`initialImb`, `isDrain`). This compares with
23 scalar fields, an integer, and three booleans in every current snapshot,
including non-WAG cases. Model-specific compact history types could avoid
carrying unused WAG slots; enabled features may permit further reduction.
The calculated Carlson/Killough fields would be restored by a pure
`recomputeDerived(config, history)` function. The current
`updateDynamicParams_()` cannot simply be called during restore: its WAG
section advances cycles and updates branch history under certain flags. The
pure recomputation function must be separate from the transition update.
It must also preserve uninitialized/sentinel behavior before the first
physical update instead of evaluating a drainage curve at a sentinel
saturation.

For `double`, the current scalar portion alone is 23 × 8 = 184 bytes per
active two-phase-law snapshot. Four scalar history values are 32 bytes;
flags and struct padding add to both figures. These are per-law snapshot
figures, not total simulator-memory savings, because the configuration and
current evaluation caches remain. Splitting out a compact non-WAG type is
what avoids allocating the unused WAG fields in each snapshot.

## Consequences of moving history out of parameters

### A problem-owned per-cell history

The smallest ownership change would put authoritative history next to the
existing explicit saturation and pressure histories in `FlowProblem`:

```text
FlowProblem
├── reservoir primary-variable time levels          existing model owner
├── per-cell accepted/candidate hysteresis history  new history owner
├── maximum saturation and minimum pressure history
└── EclMaterialLawManager
    └── curve configuration and optional derived caches
```

This would clarify the timestep transaction and make history visible to
checkpointing and diagnostics. It would also align hysteresis with
`maxOilSaturation_`, `maxWaterSaturation_`, and `minRefPressure_`.

Moving only the storage is not enough. Property evaluation currently takes a
`const MaterialLawParams&`, and the nested law reads historic values through
that parameter object. The new history must therefore reach the law through
an explicit argument or through a lightweight, per-cell evaluation view that
combines immutable configuration with const history and any calculated
coefficients. Leaving a mirrored copy of the history in parameters would
retain two authorities and most of the current synchronization burden.

For a first compatibility step, `EclMaterialLawManager` could remain the API
boundary while it reads and writes a problem-owned history store. That avoids
changing every caller at once, but the ultimate material-law evaluation API
should make both configuration and history explicit. The common material-law
code in `opm-common` should not depend on the Flow problem class itself.

### Placing history in the reservoir solution

Putting history beside the reservoir solution time levels would make acceptance
and rollback more cohesive. Putting it *inside the primary-variable vector*
would be a different and larger change. Hysteresis history is not a Newton
unknown; it should not acquire AD derivatives, matrix rows, or phase-switching
rules. The practical target is a companion history array or state bundle with
the same cell indexing and accepted/candidate lifetime, not additional PDE
unknowns.

The solution-level alternative may have cleaner ownership than `FlowProblem`
in the long run, but it touches model serialization, time-level advancement,
parallel overlap/ghost handling, and any local/distributed state transfer.
Problem ownership is a smaller first step because that class already updates
and snapshots explicit histories.

### Consequences common to either location

1. **Indexing and multiplicity.** State cannot be one scalar record per cell
   in the general case. It is per cell, active two-phase system, and possibly
   direction. The owner must preserve the same base/X/Y/Z lookup and ghost
   update rules as `EclMaterialLawManager`.
2. **Const evaluation.** Relperm and capillary-pressure calculations should
   consume a const history view at the current accepted state. Only the
   timestep-history transition mutates it. This makes it clear that a Newton
   residual evaluation does not advance hysteresis.
3. **Cache invalidation.** Changes in history still affect intensive
   quantities and their derivatives. Moving the owner does not eliminate
   invalidation; it permits a precise changed-cell or generation-based rule
   instead of the present unconditional enabled-hysteresis invalidation.
4. **Restart compatibility.** The existing internal serializer and ECL
   restart arrays have different coverage. A new state schema needs a migration
   reader or an explicit compatibility decision. WAG, PC turning points,
   directional state, and calculated coefficients need individual checks.
5. **Performance.** Recomputing inverse curve lookups and Killough coefficients
   at every relperm evaluation could be costly. A good first design keeps
   derived coefficients as a cache rebuilt once after history changes or
   restore. Snapshot reduction and removal of evaluation caches are separate
   decisions.
6. **Tests and behavior.** Historical output extrema, trapped saturation,
   relperm, capillary pressure, and Newton derivatives should match across
   timestep acceptance, a chopped timestep, in-memory restore, and restart.
   Compare full trajectories for Carlson, Killough model 4, PC hysteresis,
   WAG reversals, two-/three-phase cases, and directional settings.

## Recommended sequence

1. Introduce two phase law types or views: `HysteresisHistory` for the
   authoritative fields and `HysteresisDerived` for coefficients. Keep the
   current parameter layout initially and make `captureState()` capture only
   history. Add a pure `recomputeDerived()` called after restore.
2. Reduce non-WAG snapshots first. Test evaluation and derivative equivalence
   after every transition, not just after a single saturation excursion.
   Preserve sentinels and branch selection.
3. Treat WAG separately. Record each transition's before/after history and
   generated scanning curve. Remove a WAG field only when its reconstruction
   works for primary drainage, first imbibition, secondary drainage, later
   cycles, and rollback at each point.
4. Audit both serialization paths and output. The current `serializeOp()`
   field list is smaller than `captureState()`, and the standard restart
   arrays are smaller still. Decide which data each format promises to retain.
5. Move the resulting authoritative history to a per-cell store in
   `FlowProblem` or a companion reservoir-state object. Adapt the manager as
   a temporary bridge, then make law evaluation take explicit const history.
6. Measure memory copied per timestep and property-evaluation time before
   deciding whether calculated coefficients should remain cached.

This sequence makes snapshot size and clarity a measurable partial result.
It also produces a precise inventory for a later ownership move, without
requiring all material-law callers to change at once.

## Source map

- `opm-common/opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLawParams.hpp`:
  dynamic-state struct, history update, derived coefficient update, snapshot,
  and serializer.
- `opm-common/opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLaw.hpp`:
  values read during relperm and capillary-pressure evaluation.
- `opm-common/opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp`:
  per-cell/directional storage, update, snapshot/restore, and serialization.
- `opm-common/opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp`
  and `EclTwoPhaseMaterial.hpp`: saturation mapping, restart extrema, and
  history propagation.
- `opm-simulators/opm/simulators/flow/FlowProblem.hpp` and
  `FlowProblemBlackoil.hpp`: timestep capture/update/restore and cache
  invalidation.
- `opm-simulators/opm/simulators/flow/OutputBlackoilModule.hpp`:
  output and ECL restart mapping.
- `opm-simulators/opm/models/utils/basicparameters.hh`:
  `EnableStateRollback` default.
- `opm-common/tests/material/test_hysteresis.cpp` and
  `opm-simulators/tests/run-all-state-rollback-tests.py`:
  existing targeted rollback coverage and replay scenarios.
