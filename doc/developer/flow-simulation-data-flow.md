# OPM Flow simulation-time data flow

## Status and scope

This document is the first, reservoir-centred pass over Flow's simulation-time
data flow. It describes the default black-oil TPFA, Newton, CPU ISTL path. It
does not describe deck parsing, initialisation, grid construction, partition
construction, TPSA, NLDD, compositional Flow, GPU assembly, or reservoir
coupling except where those variants leave visible seams in the default path.

The well/group/network/gas-lift system is treated as one subsystem. Its
external contract is described here; its internal control and state
transitions are covered in the companion
[well, group, network, and gas-lift data-flow
document](flow-well-group-network-gaslift-data-flow.md).

The source snapshot inspected was:

- `opm-simulators`: `7df63fb1f821`
- `opm-common`: `8f6a6adec072`

The conclusions are based primarily on source inspection. Two serial runs with
the available `flow 2026.10-pre` binary completed successfully as smoke tests:

- `SPE1CASE2`, with 300 active cells, standard wells, and adaptive
  timestepping;
- `model2/7_HYSTERESIS_MODEL2.DATA`, with 2794 active cells and hysteresis.

A second `SPE1CASE2` run with full report steps attempted initially and the
Newton limit reduced to two iterations repeatedly exercised the convergence
failure and timestep-chop/retry path.

## Executive summary

Flow has a clear numerical pipeline, but it does not have one correspondingly
clear simulation-state object or one transaction boundary. Instead, a
conceptual timestep state is distributed across the discretisation model,
`FlowProblem`, the material-law manager, wells, aquifers, the nonlinear
system, and the linear solver.

The main data path is:

```text
current and previous primary variables
        |
        v
intensive quantities (fluid state, relperm, PVT, mobility, porosity, ...)
        |
        v
TPFA accumulation, flux, boundary, aquifer, and well source assembly
        |
        v
reservoir Jacobian and residual
        |
        v
well Schur-complement modification
        |
        v
ISTL prepare and solve
        |
        +--------------------------+
        |                          |
        v                          v
well unknown recovery       reservoir primary-variable update
        |                          |
        +-------------> next Newton iteration
```

The most important findings are:

1. **The reservoir solution is not the complete physical state.** The two
   primary-variable time levels are owned by `FvBaseDiscretization`, while
   hysteresis, irreversible rock history, VAPPARS history, mixing limits,
   well/group state, and aquifer state live in other owners.

2. **Commit and rollback are protocols, not objects.** Reservoir rollback is
   `solution(0) = solution(1)`. Well rollback is `resetWGState()`. Aquifer
   history is committed in `endTimeStep()`. Hysteresis and other explicit
   reservoir histories are advanced at the beginning of the next attempt.
   Correctness depends on lifecycle call order.

3. **`Simulator` and `Problem` act as service locators.** Components routinely
   reach through chains such as `simulator_.model().linearizer()` or
   `simulator_.problem().wellModel()`. Concrete dependencies and mutation are
   therefore much broader than method signatures suggest.

4. **Assembly and well coupling cross the same boundary several times.** Wells
   are prepared before reservoir assembly, inject connection source terms
   during TPFA assembly, contribute convergence status, eliminate well unknowns
   into the reservoir system after assembly, and reconstruct those unknowns
   after the reservoir solve.

5. **Caches participate in timestep semantics.** Intensive quantities and
   storage caches are not merely performance details: they are shifted between
   time levels, explicitly invalidated after updates, and used to represent the
   old accumulation term.

6. **The linear-solver interface exposes protocol residue.** `prepare()`
   already receives matrix and right-hand side, yet callers also invoke
   `setResidual()`, which is a no-op in `ISTLSolver`. The abstract interface
   itself contains a TODO to remove `setResidual()` and `setMatrix()`.

7. **Hysteresis is the clearest example of misplaced dynamic state.** Its
   per-cell history is nested inside material-law parameter objects, alongside
   static constitutive configuration. It needs separate capture/restore and
   restart treatment and forces property updates to invalidate intensive
   quantities.

## Concrete component and ownership graph

The compile-time property system selects the concrete types, but the runtime
ownership is approximately:

```text
Simulator
├── FIBlackOilModel
│   ├── solution_[0]                     current Newton iterate
│   ├── solution_[1]                     accepted/start-of-step state
│   ├── intensiveQuantityCache_[time]
│   ├── storageCache_[time]
│   ├── TpfaLinearizer
│   │   ├── Jacobian
│   │   └── residual
│   └── BlackOilNewtonMethod
│       └── ISTLSolverRuntimeOptionProxy / ISTLSolver
│
├── FlowProblemBlackoil
│   ├── EclMaterialLawManager
│   │   └── per-cell MaterialLawParams
│   │       └── hysteresis dynamic state
│   ├── explicit reservoir-history arrays
│   │   ├── maxOilSaturation_
│   │   ├── maxWaterSaturation_
│   │   ├── minRefPressure_
│   │   ├── polymer and mixing history
│   │   └── drift_
│   ├── BlackoilWellModel
│   │   ├── active_wgstate_
│   │   ├── last_valid_wgstate_
│   │   ├── nupcol_wgstate_
│   │   ├── per-well equation systems
│   │   ├── cellRates_
│   │   ├── network_
│   │   └── gaslift_
│   └── BlackoilAquiferModel
│       └── individual aquifer dynamic states
│
└── SimulatorFullyImplicit
    ├── AdaptiveTimeStepping
    └── NonlinearSolver
        └── NonlinearSystemBlackOilReservoir
            ├── nonlinear convergence history
            ├── relaxation history and previous update
            └── references to Simulator and BlackoilWellModel
```

This graph exposes two different meanings of “model”:

- `FIBlackOilModel` is the discretisation/state/cache owner.
- `NonlinearSystemBlackOilReservoir` is the Flow-specific nonlinear
  orchestration and convergence object.

Both are commonly reached via a method called `model()`, depending on the
calling context. This naming collision adds substantial cognitive load.

## Timestep and Newton lifecycle

### Report-step entry

`SimulatorFullyImplicit::runStep()` establishes the report step, lazily creates
one `NonlinearSolver`, invokes `beginReportStep()`, and delegates the substep
loop to `AdaptiveTimeStepping`.

The nonlinear solver owns a `NonlinearSystemBlackOilReservoir`. That object is
reused after its first creation; a nearby TUNING comment saying the solver is
re-created each report step is stale relative to the current ownership code.

### Start of a timestep attempt

`NonlinearSolver::step()` first calls
`NonlinearSystemBlackOilReservoir::prepareStep()`, whose base implementation
does the following:

```text
if previous attempt failed:
    FlowProblem::updateFailed()
        FIBlackOilModel::updateFailed()
            solution[0] = solution[1]
            invalidate/recompute current intensive quantities
else:
    FlowProblem::advanceTimeLevel()
        solution[1] = solution[0]
        shift storage cache
        shift intensive-quantity cache

set simulation time and dt
reset NewtonIterationContext
FlowProblem::beginTimeStep()
```

`FlowProblem::beginTimeStep()` then:

1. Updates explicit reservoir histories from the restored/accepted current
   state: hysteresis, maximum saturations, minimum pressure, mixing limits, and
   rock-compaction multipliers as enabled.
2. Invalidates/recomputes intensive quantities when those histories change.
3. Updates boundary-condition data when needed.
4. Calls `beginTimeStep()` on wells, aquifers, tracers, and the sequential
   temperature model.

This order is the effective state transaction. There is no object representing
“accepted timestep state”; instead, each subsystem follows its own convention.

### One Newton iteration

The standard Newton iteration is:

```text
NonlinearSystemBlackOilReservoir::nonlinearIteration()
|
├── initialLinearization()
│   ├── FlowProblem::beginIteration()
│   │   ├── BlackoilWellModel::beginIteration()
│   │   │   └── assemble(dt)
│   │   └── aquiferModel.beginIteration()
│   ├── TpfaLinearizer::linearizeDomain()
│   └── FlowProblem::endIteration()
│
├── reservoir convergence calculation
├── well/group convergence calculation
│
├── if not converged:
│   ├── BlackoilWellModel::linearize(J, r)
│   ├── ISTLSolver::prepare(J, r)
│   ├── ISTLSolver::solve(dx)
│   ├── BlackoilWellModel::postSolve(dx)
│   ├── optional nonlinear relaxation of dx
│   └── updateSolution(dx)
│       ├── BlackOilNewtonMethod::applyUpdate(solution[0], dx)
│       └── invalidate/recompute intensive quantities[0]
│
└── advance NewtonIterationContext
```

The convergence check occurs before well Schur elimination and the linear
solve. This permits an early return when the assembled reservoir and well
systems are already converged.

### Successful timestep

On convergence, `AdaptiveTimeStepping` calls `FlowProblem::endTimeStep()`.
Important state changes include:

- well postprocessing followed by `commitWGState()`;
- aquifer `endTimeStep()`, which commits cumulative influx/history;
- flow information update for output;
- drift residual capture when enabled;
- tracer and temperature postprocessing.

The reservoir primary variables are not copied into `solution[1]` here. That
happens in `advanceTimeLevel()` at the beginning of the next timestep attempt.

### Failed timestep and retry

When solve or acceptance fails, `AdaptiveTimeStepping` sets
`lastStepFailed=true`, reduces the timestep, and calls `solver.step()` again.
The next `prepareStep()` performs the rollback:

```text
reservoir: solution[0] <- solution[1]
wells:     beginTimeStep() calls resetWGState()
aquifers:  no accepted cumulative history was written because endTimeStep()
           was not called; beginTimeStep() refreshes previous connection pressure
histories: recomputed from restored current intensive quantities
caches:    current intensive quantities invalidated/rebuilt
```

This is correct only because tentative state mutations are either restored on
the next `beginTimeStep()` or deferred until `endTimeStep()`. That requirement
is not encoded in a common interface.

## Reservoir state and property evaluation

### Primary variables

`FvBaseDiscretization::solution_` is a two-slot history for implicit Euler:

- time index 0: current Newton iterate/end-of-step candidate;
- time index 1: previous accepted/start-of-step state.

Each `BlackOilPrimaryVariables` block contains numeric values plus variable
meaning, for example whether the gas switching variable represents gas
saturation, dissolved gas, or vaporised oil. Consequently, interpretation
metadata is part of the canonical solution state, not just the values.

### Primary variables to intensive quantities

For each cell, `BlackOilIntensiveQuantities::update()` constructs derived
physical state in roughly this order:

```text
primary variables and variable meanings
    -> temperature and salt
    -> phase saturations
    -> relperms from FlowProblem/material-law parameters
    -> capillary pressures and phase pressures
    -> Rs/Rv/Rsw/Rvw and saturation limits
    -> inverse formation-volume factors and viscosities
    -> phase mobilities
    -> phase densities
    -> porosity and rock-compaction multiplier
    -> extension-specific quantities
    -> flux/diffusion/dispersion quantities
```

The result owns a derived `fluidState_`, mobilities, porosity, and module
quantities. It is cached by cell and time index in the discretisation model.

The evaluation's apparent inputs are misleadingly small. Through `Problem` it
also reads:

- material-law parameters and saturation-region information;
- static and dynamic rock properties;
- maximum historic saturations and minimum historic pressure;
- composition-change limits;
- schedule/report-step state;
- linearisation mode;
- module configuration;
- the global simulator in some extension paths.

The cache must be invalidated after every primary-variable update and after a
dynamic property-history update.

### Static, dynamic, derived, and workspace data

The present type structure does not consistently separate four categories:

1. Static constitutive configuration, such as table selection and endpoint
   scaling configuration.
2. Accepted dynamic physical history, such as hysteresis turning points,
   maximum saturation, and irreversible minimum pressure.
3. Derived state, such as fluid properties and mobilities.
4. Assembly workspace, such as AD evaluations and cached storage blocks.

That makes it difficult to know whether passing a “parameter” object is
read-only, whether a cached object is safe to reuse, and which pieces must be
rolled back or serialised.

## Hysteresis case study

Hysteresis is enabled through the material-law type and configured by
`EclHysteresisConfig`, but the dynamic per-cell state is stored inside
`EclHysteresisTwoPhaseLawParams`. Examples include historic extrema, scanning
curve quantities, and WAG process state.

The simulation lifecycle is:

```text
accepted/restored solution[0]
    -> cached current intensive quantities
    -> FlowProblemBlackoil::updateExplicitQuantities_()
    -> FlowProblem::updateHysteresis_()
    -> EclMaterialLawManager::updateHysteresis(cell params, fluid state)
    -> mutation of MaterialLawParams
    -> invalidate/recompute current intensive quantities
    -> subsequent relperm/capillary-pressure evaluation reads mutated params
```

Consequences:

- A nominal parameter manager owns accepted physical history.
- The same per-cell object mixes static curves/configuration and mutable state.
- Updating history changes future property evaluations, so it invalidates a
  cache owned by another object.
- Restart serialization must specially serialize dynamic contents of every
  material-law parameter object.
- Output extracts hysteresis state through the material-law manager instead of
  reading the reservoir solution.
- `EclMaterialLawManager` has explicit begin-timestep capture/restore support,
  but no production call site currently uses it. Flow instead relies on
  hysteresis being updated only from accepted/restored state at timestep
  entry.
- Directional material-law parameters multiply the capture/restore surface.

A simpler target is:

```text
ReservoirHistoryState[cell]
├── hysteresis state
├── maximum saturations
├── minimum irreversible pressure
└── other accepted constitutive history

MaterialLawConfiguration[cell or region]
└── immutable curves, endpoints, options, and region mappings
```

Property evaluation would then accept a const configuration and a const state
view. Committing or rolling back a timestep would move/copy/swap one explicit
history state together with the primary variables.

## TPFA assembly data flow

`TpfaLinearizer` owns the global Jacobian and residual. Each global assembly
resets them, then performs:

1. **Interior face fluxes.** Current intensive quantities from both cells and
   transmissibility/face data produce an AD flux. Values enter the residual;
   derivatives enter diagonal and neighbour Jacobian blocks.
2. **Accumulation.** Current storage is computed from time-index-0 intensive
   quantities. Previous storage comes from time-index-1 intensive quantities
   or from the shifted storage cache. Their difference is scaled by volume/dt.
3. **Boundary contributions.** Boundary fluxes add residual and Jacobian
   blocks.
4. **Dense cell sources.** Ordinary source terms and aquifer sources are
   evaluated through `FlowProblem::source()`.
5. **Sparse sources.** In the normal separated-source path,
   `BlackoilWellModel::addReservoirSourceTerms()` inserts connection rates and
   their local derivatives into reservoir rows.

The assembly boundary is not pure:

- it reads model caches and `FlowProblem` dynamic history;
- aquifer source evaluation updates per-connection current pressure/rate
  workspace;
- well source evaluation reads a well equation system assembled immediately
  before reservoir linearisation;
- storage caches may be updated during first-iteration assembly;
- flow and velocity output buffers may also be updated.

### Well elimination after reservoir assembly

After convergence has been checked and a solve is required,
`BlackoilWellModel::linearize(J, r)` applies the well connection auxiliary
module. Per well it:

- optionally adds explicit matrix connectivity/contributions;
- gathers residual blocks for perforated cells;
- applies the eliminated well system to the residual;
- scatters the modified blocks back.

For a standard well this is the Schur operation represented by the stored
well blocks, conceptually

```text
J_res <- J_res - C^T D^-1 B
r_res <- r_res - C^T D^-1 r_well
```

After solving the reduced reservoir system, `postSolve(dx)` gathers the
reservoir updates at each well's perforation cells, recovers well unknowns,
and updates active well state.

## Well/group/network/gas-lift external contract

The first-pass boundary of the subsystem is:

| Phase | Reservoir-to-well input | Well-to-reservoir/output effect |
|---|---|---|
| Report-step entry | schedule, summary state, report-step index | well/group/network structures |
| Timestep entry | accepted cell pressures/properties, dt | reset active state from last valid state; initialise wells and controls |
| Iteration entry | current reservoir intensive quantities, iteration context | assembled well equations, controls/network update, connection rates |
| Reservoir assembly | none beyond already assembled well data | connection source residual/Jacobian terms |
| Convergence | reservoir `B_avg` and reservoir convergence status | well, group-control, and network failures |
| Linearisation | reservoir Jacobian and residual | Schur-complement modification |
| Post solve | reservoir update at perforation cells | recovered well update; mutation of active well/group state |
| Successful step | accepted reservoir state, time, dt | potentials, limits, cumulative state; commit active to last-valid state |
| Failed retry | retry lifecycle only | next `beginTimeStep()` resets active from last-valid state |

The external interface already shows why this subsystem needs a separate deep
investigation: “well state” is not one object or one phase. The generic model
stores active, last-valid, and NUPCOL `WGState` values; the concrete model also
owns well-local linear systems, cell-rate caches, network and gas-lift helpers,
and temporary reconstructed solutions.

## Linear solver boundary

The default call path is:

```text
TpfaLinearizer-owned matrix/residual
    -> well Schur modification
    -> ISTLSolver::prepare(matrix, residual)
       - retain non-owning pointers
       - modify overlap rows for some parallel configurations
       - create or update preconditioner/solver hierarchy
       - potentially construct a well-aware linear operator
    -> ISTLSolver::setResidual(residual)   [currently a no-op]
    -> ISTLSolver::solve(update)
    -> iteration count and convergence diagnostics
```

`ISTLSolver` is not independent of Flow state. Besides matrix and RHS it reads:

- grid overlap/interior information;
- simulator and model for true-IMPES CPR weights;
- iteration context for preconditioner-reuse policy;
- well model and number of active local wells;
- runtime solver parameters and MPI communication.

It may also mutate the input matrix's overlap rows, retains raw pointers to the
matrix and RHS between `prepare()` and `solve()`, and caches solver hierarchy
state across calls.

A clearer interface would combine the call protocol:

```text
LinearSystemView
    matrix
    rhs
    parallel_layout
    block_layout
    optional well_operator
    state/version identifiers

LinearSolveResult
    update
    converged
    iterations
    reduction
    setup_reused
```

The solver may still cache a preconditioner, but reuse decisions should be
based on explicit system/topology/version metadata rather than reaching back
into `Simulator`, `Problem`, and `WellModel`.

## State ownership inventory

| State | Current owner | Kind/lifetime | Main writers | Commit/rollback and serialization |
|---|---|---|---|---|
| `solution_[0]` | `FvBaseDiscretization` | canonical tentative reservoir state; Newton/timestep | Newton update, restart load | accepted later by copy to slot 1; reset from slot 1; serialized |
| `solution_[1]` | `FvBaseDiscretization` | accepted/start-of-step reservoir state | `advanceTimeLevel()` | rollback source; serialized |
| primary-variable meaning | each `BlackOilPrimaryVariables` | canonical phase-presence/switch state | Newton update switching logic | follows solution vectors |
| intensive-quantity cache | `FvBaseDiscretization` | derived per-cell/time cache | cache update/invalidation | shifted at time-level advance; reconstructed after restart |
| storage cache | `FvBaseDiscretization` | derived accumulation cache | first assembly, rebuild | shifted at time-level advance; reconstructed after restart |
| Jacobian and residual | `TpfaLinearizer` | per-Newton workspace | reservoir assembly, wells | reset each assembly; not physical restart state |
| explicit saturation/pressure history | `FlowGenericProblem` | accepted physical history | timestep-entry update functions | relies on update ordering; serialized selectively |
| hysteresis dynamic state | per-cell material-law params in `EclMaterialLawManager` | accepted physical history embedded in parameters | `updateHysteresis_()` | special parameter serialization; implicit retry semantics |
| mixing-rate controls/history | `FlowProblemBlackoil::mixControls_` and related objects | accepted physical/control history | timestep entry and schedule changes | serialized through problem-specific path |
| drift compensation | `FlowProblem::drift_` | accepted derived/history state | successful `endTimeStep()` | serialized |
| active well/group state | `BlackoilWellModelGeneric::active_wgstate_` | tentative Newton/timestep state | well assembly, controls, `postSolve()` | committed to last-valid on success; serialized |
| last-valid well/group state | `BlackoilWellModelGeneric::last_valid_wgstate_` | accepted rollback state | `commitWGState()` | copied to active on timestep entry; serialized |
| NUPCOL well/group state | `BlackoilWellModelGeneric::nupcol_wgstate_` | frozen iteration-policy state | `updateNupcolWGState()` | neither the canonical active nor accepted state |
| per-well equation blocks | individual well objects | per-Newton assembly workspace | well assembly | reused by source, Schur, and recovery phases |
| connection cell-rate cache | `BlackoilWellModel::cellRates_` | per-Newton derived cache | well assembly | rebuilt after well equations |
| well reconstructed solution cache | `BlackoilWellModel` | solve/post-solve workspace | system solve path | consumed during post-solve |
| aquifer cumulative influx/history | individual aquifer | accepted physical history | successful `endTimeStep()` | no commit on failed attempt; serialized |
| aquifer current rate/pressure | individual aquifer | iteration workspace | source evaluation | recomputed; some parts serialized for restart/output |
| residual norms and previous Newton update | `NonlinearSystem` | nonlinear-algorithm history per timestep | convergence/update stabilization | cleared at timestep init; not physical state |
| iteration counters/init flag | `FlowProblem`'s `NewtonIterationContext` | lifecycle control | nonlinear system | reset per attempt; advanced per Newton iteration |
| timestep proposal/restart count | `AdaptiveTimeStepping` and timers | timestep-control state | time-step controller | advanced on success, reduced on failure |
| TUNING/TUNINGDP model parameters | `SimulatorFullyImplicit` and `NonlinearSystem` | duplicated configuration | report-step callback | manually updated in both owners |
| preconditioner/solver hierarchy | `ISTLSolver` | cross-iteration performance state | `prepare()` | reused according to implicit simulator/well context |

The serialization column is itself diagnostic: restart state is assembled by
serializing the model, problem, material-law manager, wells, aquifers, and
extension containers separately. There is no single authoritative schema for
simulation state.

## Interface-edge inventory

| Producer -> consumer | Data crossing boundary | Mutation/hidden dependency |
|---|---|---|
| Adaptive timestepping -> nonlinear solver | time, dt, report/substep indices, failure flag | solver reads problem iteration state indirectly |
| Nonlinear system -> `FlowProblem` | lifecycle calls | calls fan out to wells, aquifers, caches, histories |
| solution -> intensive quantities | primary-variable block and time index | evaluation reaches through `Problem` to many dynamic inputs |
| material-law manager -> intensive quantities | per-cell law params, relperm, capillary pressure | nominal params contain mutable history |
| intensive quantities -> TPFA assembly | fluid state, mobility, porosity, module quantities | cache validity external to type |
| old-state cache -> accumulation | cached storage or time-index-1 IQ | first-iteration assembly may update cache |
| aquifer -> source assembly | water/energy rate as AD value | source evaluation updates aquifer workspace |
| wells -> TPFA source assembly | connection rates and derivatives | assumes well equations assembled earlier in same iteration |
| reservoir assembly -> convergence | residual, pore volumes, fluid properties | includes MPI reductions and well status |
| wells -> assembled linear system | Schur contribution | mutates matrix and RHS after reservoir assembly |
| linearizer -> ISTL | matrix and RHS by reference | solver retains pointers and may mutate overlap rows |
| ISTL -> nonlinear system | reservoir update and diagnostics | solver also consults well/model/iteration state indirectly |
| reservoir update -> wells | perforation-cell update blocks | reconstructs well unknowns and mutates active state |
| reservoir update -> model | full update vector | mutates current primary variables and phase meanings |
| accepted solution -> output/restart | solution plus distributed side-state | consumers gather from several owners |

## Why the current flow is difficult to understand

### 1. State is classified by implementation location, not physical meaning

The primary solution, constitutive history, well state, aquifer history, and
algorithm state all use different ownership and commit conventions. A reader
cannot start from a “simulation state” type and discover the whole state.

### 2. Method signatures under-report dependencies

`BlackOilIntensiveQuantities::update(problem, priVars, cell, time)` and
`ISTLSolver::prepare(matrix, rhs)` look relatively narrow, but both can access
large portions of the simulator object graph. This impedes isolated tests and
makes apparently local changes globally significant.

### 3. Temporal coupling substitutes for explicit data contracts

Examples:

- well equations must be assembled in `beginIteration()` before TPFA asks for
  connection rates;
- well source terms must be in the reservoir residual before convergence;
- well Schur elimination must occur only after convergence says a solve is
  needed;
- `postSolve()` must precede reservoir update;
- hysteresis must update only after current solution has been accepted or
  restored;
- well rollback occurs on the next attempt's `beginTimeStep()`, not at the
  failure point.

### 4. “Parameters”, “model”, and “state” are overloaded words

Material-law parameters contain dynamic state. Two distinct objects are called
models. `WellState` is one part of a wider `WGState`. Linear-solver parameters
coexist with cached solver state and simulator-derived policy.

### 5. Cache management is externally coordinated

Callers need to know when a property-history mutation invalidates intensive
quantities and when time-level advance should shift storage caches. Cache
versioning is represented by mutable validity vectors and call order rather
than by explicit dependency versions.

### 6. Compile-time polymorphism hides the runtime graph

The property system is useful for assembling variants, but a reader must
resolve `Model`, `Problem`, `Linearizer`, `LocalResidual`, `WellModel`, and
`LinearSolverBackend` before ordinary calls become concrete. Runtime control
then adds virtual auxiliary modules and solver proxies.

### 7. Some interface comments no longer match implementation

Examples include the TUNING comment about recreating the nonlinear solver each
report step and the abstract ISTL protocol requiring calls that the concrete
solver ignores. Stale protocol documentation is particularly costly in an
architecture already dependent on call ordering.

## Recommended simplification sequence

The following order creates useful seams without requiring a full rewrite.

### A. Introduce explicit state taxonomy and views

Define names and lightweight views for:

- `ReservoirPrimaryState`;
- `ReservoirHistoryState`;
- `ReservoirDerivedState`;
- `ReservoirAssemblyWorkspace`;
- `WellGroupState`;
- `AquiferState`;
- `NonlinearAlgorithmState`.

The first implementation can wrap existing storage rather than move it. The
immediate benefit is explicit method contracts and a common inventory.

### B. Make timestep transaction operations explicit

Replace distributed lifecycle assumptions with a small protocol:

```text
beginAttempt(accepted_state, dt)
commitAttempt(candidate_state)
rollbackAttempt()
```

Initially these methods can delegate to existing `advanceTimeLevel()`,
`updateFailed()`, `resetWGState()`, and subsystem hooks. The key improvement is
one auditable transaction boundary.

### C. Separate material-law configuration from dynamic history

Move hysteresis dynamic fields out of material-law parameter objects into
`ReservoirHistoryState`. Convert saturation-property evaluation to take:

```text
evaluate(configuration, history_view, fluid_state)
```

This is the highest-value ownership correction because it simplifies
rollback, restart, output, testing, and cache invalidation simultaneously.

### D. Encapsulate property evaluation behind a const cell-state interface

Create a `CellPropertyContext` that explicitly contains the current primary
state, accepted history, region indices, immutable tables, time context, and
linearisation mode. Remove `Simulator` access from the common property path.

Extension modules can migrate incrementally; the current elemCtx-less overload
already demonstrates a partial direction, although it supports only a subset
of features.

### E. Make assembly stages explicit

Represent assembly as named stages with owned inputs/outputs:

```text
prepareCoupledSources()
assembleReservoir()
evaluateConvergence()
eliminateCoupledUnknowns()
solveReducedSystem()
recoverCoupledUnknowns()
applyStateUpdate()
```

This preserves the current algorithm while making the well contract and
ordering visible.

### F. Collapse the linear-solver protocol

Replace `prepare()` + `setResidual()` + implicit simulator queries with one
`solve(LinearSystemView, ReuseContext)` operation. Return diagnostics rather
than exposing solver state through later getters.

### G. Add cache dependency versions

Associate intensive/storage cache entries with explicit versions of primary
state, constitutive history, and configuration. This would replace broad
“invalidate everything after mutation” knowledge with checkable dependencies.

### H. Rename the two model layers

Even without structural changes, naming the objects according to role would
help:

- `FIBlackOilModel` -> discretisation or reservoir-state model;
- `NonlinearSystemBlackOilReservoir` -> nonlinear reservoir driver/system.

## Candidate low-risk changes

These can be considered before moving state ownership:

1. Remove the redundant `setResidual()`/`setMatrix()` protocol from the common
   ISTL interface where implementations permit it.
2. Pass an explicit `LinearSystemView` from nonlinear system to solver.
3. Replace repeated `simulator_.model().linearizer()` chains in
   `NonlinearSystemBlackOilReservoir` with constructor-injected narrow
   references or accessors.
4. Document or type the sign convention for residual, right-hand side, and
   Newton update.
5. Put the accepted/rejected attempt transition in one named orchestration
   method and add lifecycle tests.
6. Add debug assertions for well source cache freshness and equation-assembly
   generation.
7. Correct stale solver-lifetime and linear-solver-protocol comments.
8. Add a state inventory test that round-trips serialization after an accepted
   timestep and compares all accepted-state owners.

## Higher-risk architectural changes

1. Move hysteresis and other irreversible per-cell histories next to the
   solution vectors.
2. Replace the global `Simulator` service-locator pattern in property,
   assembly, well, and solver code with explicit context objects.
3. Give reservoir, well, and aquifer candidates a common timestep transaction.
4. Separate well control/network iteration state from the well equation state
   and from accepted cumulative output state.
5. Make the coupled nonlinear system own an explicit block-system description
   rather than modifying a reservoir-only system in several phases.

## Runtime validation still worth adding

The smoke run validated the selected executable path, not every state
transition. A small opt-in lifecycle tracer would make the source conclusions
regression-testable. It should log only:

- report step, substep, attempt, and Newton iteration;
- state generation identifiers for solution, history, IQ cache, well state,
  and aquifer state;
- begin/commit/rollback events;
- assembly generation for well equations, reservoir system, and Schur update;
- matrix/RHS versions passed to the solver.

The smoke runs exercised all three scenarios below at executable level, but
they did not expose internal state generations. A permanent tracer-based
regression suite should use:

1. `SPE1CASE2` for the normal reservoir/well path.
2. `model2/7_HYSTERESIS_MODEL2.DATA` for constitutive history.
3. A deliberately forced timestep rejection for rollback.

The tracer should verify invariants rather than dump full field values. For
example, after retry entry the current reservoir generation and active well
state must descend from the same accepted attempt.

## Deferred well/group/network/gas-lift investigation

The next investigation should begin from the external contract above and
answer these questions:

1. Which fields in `active_wgstate_`, `last_valid_wgstate_`, and
   `nupcol_wgstate_` are authoritative, derived, cumulative, or temporary?
2. What exact operations mutate state during `updateWellControlsAndNetwork()`?
3. How are group targets, guide rates, network node pressures, well controls,
   and gas-lift decisions ordered and iterated?
4. Which state changes are Newton-iteration-local, timestep-tentative,
   report-step-scoped, or accepted cumulative state?
5. Which convergence failures request another Newton iteration versus an
   internal network/control iteration?
6. What state is communicated across ranks, and who owns the authoritative
   value for distributed wells and groups?
7. Can the subsystem expose a single coupled-equation contract containing
   source terms, convergence, Schur elimination, recovery, commit, and
   rollback?
8. Can guide-rate, group-control, network, and gas-lift state be separated from
   the well equation objects without changing numerical behaviour?

That work should produce its own internal state inventory and sequence
diagrams. It should not be appended to the reservoir diagram until its
transaction boundaries and internal iteration nesting are explicit.

## Primary source map

The main files followed in this investigation are:

- `opm/simulators/flow/SimulatorFullyImplicit_impl.hpp`
- `opm/simulators/timestepping/AdaptiveTimeStepping_impl.hpp`
- `opm/simulators/flow/NonlinearSolver.hpp`
- `opm/simulators/flow/NonlinearSystem_impl.hpp`
- `opm/simulators/flow/NonlinearSystemBlackOilReservoir_impl.hpp`
- `opm/simulators/flow/FlowProblem.hpp`
- `opm/simulators/flow/FlowProblemBlackoil.hpp`
- `opm/models/discretization/common/fvbasediscretization.hh`
- `opm/models/discretization/common/tpfalinearizer.hh`
- `opm/models/blackoil/blackoilprimaryvariables.hh`
- `opm/models/blackoil/blackoilintensivequantities.hh`
- `opm/models/blackoil/blackoillocalresidualtpfa.hh`
- `opm/simulators/wells/BlackoilWellModel.hpp`
- `opm/simulators/wells/BlackoilWellModel_impl.hpp`
- `opm/simulators/wells/BlackoilWellModelGeneric.hpp`
- `opm/simulators/wells/WellConnectionAuxiliaryModule.hpp`
- `opm/simulators/linalg/AbstractISTLSolver.hpp`
- `opm/simulators/linalg/ISTLSolver.hpp`
- `opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp`
- `opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLawParams.hpp`
