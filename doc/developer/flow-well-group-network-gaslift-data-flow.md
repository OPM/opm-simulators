# OPM Flow well, group, network, and gas-lift data flow

## Status and scope

This document is the second pass over Flow's simulation-time data flow. It is
a companion to [flow-simulation-data-flow.md](flow-simulation-data-flow.md)
and examines the internals of the well, group-control, production-network, and
gas-lift subsystem. Particular attention is paid to the points where the
reservoir simulator enters the subsystem and where data returns to reservoir
assembly, convergence, the linear solve, restart, and output.

The main path considered is black-oil Flow with TPFA, Newton, standard or
multisegment wells, and the extended production network. This is a structural
investigation rather than a description of every supported control keyword or
the detailed mathematics of every well and gas-lift algorithm. Reservoir
coupling, NLDD, GPU-specific extraction, polymer/filter-cake details, well
testing, and actions are mentioned only where they expose an important
interface or state-management property.

The source snapshot inspected was:

- `opm-simulators`: `e015dec306f6`
- `opm-common`: `c68d905fcdc4`

The conclusions are based on source inspection. No new simulation result is
used as evidence in this pass.

## Executive summary

The subsystem has one obvious facade, `BlackoilWellModel`, but it does not
have one narrow input/output contract. During one Newton iteration it reads
the reservoir through `Simulator`, changes well and group state, repeatedly
solves local well equations, balances the network, possibly optimizes gas
lift, contributes sources and a Schur complement to the reservoir equations,
and finally reconstructs well unknowns from the reservoir update.

The central feedback cycle is:

```text
WellState rates/potentials and GroupState controls
          |
          v
aggregate group rates and allocate group targets
          |
          v
switch well/group controls
          |
          v
compute network node pressures --> impose well THP limits
          |                              |
          |                              v
          |                     solve local well equations
          |                              |
          +------------------------------+
          |
          v
gas-lift optimization --> change ALQ and well potentials
          |
          v
refresh guide rates and, when required, repeat the cycle
```

The main findings are:

1. **`BlackoilWellModel` is both a component and an orchestrator.** It owns
   state, constructs timestep-local well objects, coordinates group control,
   network, and gas lift, and implements the numerical reservoir/well
   coupling. This gives its public surface several unrelated roles.

2. **The important boundary is crossed multiple times per Newton
   iteration.** The reservoir supplies fluid properties through the
   `Simulator`; wells supply connection source terms during TPFA assembly;
   well residuals affect convergence; eliminated well equations modify the
   reservoir matrix and residual; and the reservoir update is used to recover
   the well update.

3. **The facility algorithms communicate primarily by shared mutable
   state.** `WellState` and `GroupState` are passed widely and are modified by
   target allocation, control switching, network constraints, local well
   solves, and gas-lift optimization. Return values usually report only that
   something changed, not what changed.

4. **There are three copies of well/group/test state with different temporal
   meanings.** `active_wgstate_` is tentative, `last_valid_wgstate_` is the
   accepted rollback point, and `nupcol_wgstate_` freezes selected values after
   the NUPCOL window. Network pressures have their own active/accepted pair.
   Correctness depends on knowing which copy each calculation reads.

5. **The production network is an explicit fixed-point iteration outside the
   global Jacobian.** Rates determine node pressures, node pressures become
   well THP limits, local well solves change rates, and the cycle repeats.
   Damping and an outer/inner iteration protocol provide convergence.

6. **Gas lift is a controller inside that network fixed point.** It consumes
   network pressures, well potentials, group constraints, and guide-rate
   information, then mutates each well's ALQ and potentials. An ALQ change is
   a reason to rebalance the network and update guide rates.

7. **Group data is both state and derived cache.** `GroupState` contains
   control modes and history, but also globally aggregated rates, reductions,
   targets, counts, and network leaf rates. Much of it is recomputed and
   communicated repeatedly, yet it shares one type and serialization protocol
   with longer-lived control state.

8. **Parallel ownership leaks into the domain model.** A distributed well can
   exist on several ranks, while `WellState` also carries global information
   for all wells. Well equation blocks, mixing rates, group rates, control
   changes, network imbalance, and gas-lift allocations each have different
   collective communication rules.

9. **Rollback is mostly coherent but still protocol-driven.** The accepted
   `WGState` and network pressure state are restored at the next timestep
   attempt. However, timestep-local well objects also contain mutable
   operability, dynamic limits, equation caches, and switching logs that are
   rebuilt or reinitialized by a separate sequence. The pre-step network
   rebalance itself contains a TODO questioning backup behavior.

10. **A simpler interface should begin with explicit state snapshots and
    phase results, not by merging more classes.** Separating accepted physical
    state, Newton-iterate state, NUPCOL policy state, and derived aggregates
    would make the existing algorithms easier to reason about before changing
    them.

## Ownership and state graph

The principal runtime ownership is approximately:

```text
FlowProblemBlackoil
└── BlackoilWellModel
    ├── active_wgstate_                         tentative timestep state
    │   ├── WellState
    │   │   └── SingleWellState per local well
    │   │       ├── bhp, thp, rates, potentials, control modes
    │   │       ├── perforation and segment state
    │   │       ├── group targets and retained network THP limit
    │   │       └── ALQState
    │   ├── GroupState
    │   │   ├── group rates and reductions
    │   │   ├── production/injection control modes
    │   │   ├── REIN, VREP, sales, and pressure-maintenance data
    │   │   └── network leaf rates and auto-choke THP
    │   └── WellTestState
    ├── last_valid_wgstate_                     accepted/rollback copy
    ├── nupcol_wgstate_                         iteration-policy snapshot
    ├── well_container_                         timestep-local WellInterface objects
    │   ├── StandardWell or MultisegmentWell
    │   ├── primary variables and operability state
    │   └── local equation blocks B, C, D and well residual
    ├── GroupStateHelper                        broad service facade
    ├── GuideRate and GuideRateHandler
    ├── BlackoilWellModelNetwork
    │   ├── current and accepted node pressures
    │   ├── current and accepted branch data
    │   └── auto-choke group-THP cache
    ├── BlackoilWellModelGasLift
    │   └── last optimization time
    ├── VFPProperties and RateConverter
    └── other well histories and caches
```

The schedule, summary state, VFP tables, reservoir intensive quantities, and
parallel communicator are not owned by this graph. They are reached through
references to `Simulator`, `Schedule`, or `BlackoilWellModelGeneric`.

There are also two representations of an active well:

- `SingleWellState` is the persistent, serializable simulation state.
- `StandardWell` or `MultisegmentWell` is a timestep-local computational
  object containing equations, primary-variable representation, control and
  operability machinery, and references to schedule/static data.

Keeping persistent state outside the polymorphic numerical object is useful
for reconstruction and distributed wells, but the split makes it necessary to
know whether a value lives in `SingleWellState`, in `WellInterface`, or in
both. Dynamic network THP is a concrete example: it is applied to the well
object and retained in `SingleWellState::network_thp_limit` so it survives
well reconstruction or network detachment.

The important state stores and lifetimes are:

| Store | Meaning and lifetime | Commit/restore behavior |
| --- | --- | --- |
| `active_wgstate_` | Tentative well, group, and well-test state for the current timestep attempt | Mutated throughout facility calculations; copied from `last_valid_wgstate_` on attempt entry |
| `last_valid_wgstate_` | Last accepted well/group/test state | Replaced at report-step setup, selected topology changes, and successful timestep completion |
| `nupcol_wgstate_` | Snapshot used to freeze selected control-allocation inputs | Refreshed during the NUPCOL window and for specified REIN/VREP rate changes; serialized separately |
| `well_container_` | Timestep-local numerical well objects | Recreated at timestep entry; not the primary restart/rollback representation |
| Per-well B/C/D and residual | Linearization at the latest well assembly point | Cleared/reassembled; valid only for the matching reservoir/control state |
| `cellRates_` | Per-cell view of latest connection rates | Rebuilt after well equation assembly |
| `guideRate_` | Durable guide-rate values and timing policy | Serialized, but not part of WG commit/reset; normally refreshed during timestep entry and when invalidated |
| Network current maps | Tentative node pressures and branch data | Iterated and damped during facility updates |
| Network last-valid maps | Accepted node pressures and branch data | Committed/restored alongside WG state by separate calls |
| `well_group_thp_calc_` | Auto-choke root-finding cache | Serialized with the network; not included in network commit/reset copying |
| `last_glift_opt_time_` | Gas-lift scheduling history | Serialized, but has no accepted/tentative pair; retry at the same simulation time is explicitly eligible to optimize again |
| VFP/rate-conversion objects | Report-step computational services derived from schedule/reservoir state | Rebuilt/redefined rather than copied as timestep state |

## External entry and exit points

The table below is the short version of the subsystem's contract with the rest
of Flow.

| Simulation phase | Entry into subsystem | Data entering | State changed internally | Data returning |
| --- | --- | --- | --- | --- |
| Report-step start | `beginReportStep()` | Schedule wells/groups, reservoir state, restart state | Local well structure, `WellState`, group structure, VFP and rate-conversion objects; accepted snapshot | None directly |
| Timestep-attempt start | `beginTimeStep()` | Accepted WG/network state, schedule events, reservoir properties, time and dt | Restores state, creates well objects, potentials, guide rates, targets, optional initial well solves | Prepared well/group/network state |
| Newton begin | `beginIteration()` -> `assemble(dt)` | Current reservoir intensive quantities and iteration context | Controls, targets, network pressure, ALQ, well primary variables and equation systems | Cell connection-rate cache and change flags |
| Reservoir TPFA assembly | `addReservoirSourceTerms()` / `computeTotalRatesForDof()` | Reservoir residual/Jacobian storage or cell index | Normally only assembly caches | Perforation source residuals and diagonal derivatives |
| Convergence | `getWellConvergence()` | Average formation-volume factors and reservoir-converged flag | Convergence report/logging | Well failures, control-target violation, forced network iteration |
| Before linear solve | `linearize(J, r)` | Assembled reservoir Jacobian and residual | Well elimination scratch | Schur-complement matrix and residual modification |
| Optional preconditioner/GPU paths | `addWellPressureEquations()`, `addBCDMatrix()`, `getWellContributions()` | Reservoir pressure weights or extraction buffers | None conceptually | Well blocks in solver-specific forms |
| After linear solve | `postSolve(dx)` | Reservoir Newton update at perforated cells | Recovers well unknown update and mutates `WellState` | Updated well state for next Newton iteration |
| Successful timestep | `timeStepSucceeded()` | Converged reservoir/well state, time and dt | Potentials, PI, test/economic state, histories; commits WG/network state | State for output and next timestep |
| Output | `wellData()`, `groupAndNetworkData()`, WBP methods | Current state and schedule metadata | Mainly report-object construction; network also computes an undamped reporting solution | `data::Wells`, group controls/guide rates, node and branch data |
| Restart serialization | `FlowProblem::serializeOp()` -> well-model serialization | Serializer | None | Active, accepted, and NUPCOL WG states; network and gas-lift state; other well histories |

Two points deserve emphasis:

- Reservoir assembly does not pass a compact `ReservoirView` into the well
  subsystem. Individual well algorithms reach reservoir intensive quantities,
  time, iteration context, schedule, summary state, and linearizer services
  through `Simulator` and `BlackoilWellModel`.
- The well subsystem does not return one `WellAssemblyResult`. Its outputs are
  delivered through source-term callbacks, mutation of the global Jacobian and
  residual, convergence flags, and persistent state mutation at different
  times.

## Lifecycle in detail

### Report-step entry

`BlackoilWellModel::beginReportStep()` performs the schedule/topology-level
transition:

1. It records the report step and well/group events.
2. It constructs a new `RateConverter`.
3. It initializes the rank-local scheduled-well structure and parallel-well
   information. Unless a deserialized state is pending, it also initializes
   perforation and well state.
4. It initializes group control modes and any GPMAINT regional-pressure
   calculators.
5. It defines rate-conversion state from the reservoir and rebuilds VFP
   properties.
6. It commits the resulting `WGState` as a rollback point.

This is also where the persistent/computational split begins: schedule-derived
well structure and `WellState` are initialized here, but `well_container_` is
created at timestep start.

### Start or retry of a timestep

`BlackoilWellModel::beginTimeStep()` first handles a possible topology change
caused within a report step, then calls `resetWGState()`. The reset copies
`last_valid_wgstate_` into `active_wgstate_`, restores accepted network node
and branch data, and retargets `GroupStateHelper` to the restored active
objects.

It then performs a broad preparation sequence:

```text
restore accepted WG/network state
  -> initialize default ALQ and gas-lift timestep counters
  -> test wells
  -> create StandardWell/MultisegmentWell objects
  -> aggregate initial group data
  -> initialize each well object and network-imposed limits
  -> attach VFP and guide-rate services
  -> initialize well-specific histories and completion status
  -> compute event-driven well potentials
  -> update guide rates
  -> update pressure-maintenance targets
  -> aggregate group data and allocate group targets
  -> solve new or materially changed wells for an initial state
```

The sequence shows that `beginTimeStep()` is not just reset/initialization. It
already runs domain algorithms and can mutate rates, controls, primary
variables, potentials, and group data before the first global Newton
iteration.

### One Newton iteration: facility preparation

The reservoir calls `BlackoilWellModel::assemble(dt)` from
`FlowProblem::beginIteration()`. On the first required initialization within a
timestep, it calculates explicit well quantities and performs timestep
preparation. The central call is then `updateWellControlsAndNetwork()`.

Its outer fixed-point iteration calls
`updateWellControlsAndNetworkIteration()`:

```text
1. updateAndCommunicateGroupData(update targets = true)
2. updateWellControls()
3. network.update()
4. gaslift.maybeDoGasLiftOptimize()
5. prepareWellsBeforeAssembling()
6. if required, update guide rates
7. decide whether network/ALQ feedback requires another outer iteration
```

After the outer loop, `assemble()` performs one non-iterating assembly of the
well equations and creates `cellRates_` from connection rates for reservoir
source lookup.

The order is significant. Group targets are calculated from the rates and
guide-rate state available at the start of an outer iteration. Network
pressures then modify well THP constraints. Gas lift can modify ALQ and
potentials. Local well solves update rates under those choices. If ALQ or the
network remains unsettled, the next outer iteration observes the new rates and
recomputes group data and targets.

`well_group_control_changed` is not itself part of the condition that repeats
this network outer loop. It is carried into the convergence report so the
global Newton process performs another iteration. Network imbalance and ALQ
changes determine the local outer-loop continuation.

If the configured maximum network outer iterations is reached while the
network is expected to balance on the next Newton iteration,
`network_needs_more_balancing_force_another_newton_iteration_` is set. Thus
network fixed-point convergence can deliberately control global Newton
convergence without contributing a network equation to the global system.

### Local well solve versus global well coupling

`prepareWellsBeforeAssembling()` calls
`WellInterface::prepareWellBeforeAssembling()` for each local well. During the
early Newton iterations this commonly runs an inner nonlinear well solve:

```text
copy/update well primary variables
  -> assemble local well residual and B/C/D derivative blocks
  -> test local convergence
  -> solve D * delta_w = residual_w
  -> update SingleWellState and primary variables
  -> optionally switch status or control
  -> repeat
```

This is an approximate/local solve with reservoir unknowns held at their
current values. It improves the initial state and handles controls and
operability, but it does not replace the global reservoir/well coupling.

After all control/network/gas-lift work, the final
`assembleWellEqWithoutIteration()` builds the equation blocks at the state
that will be coupled to the reservoir. For a standard well the conceptual
linearized system is:

```text
[ A  C^T ] [delta_r] = [r_r]
[ B   D  ] [delta_w]   [r_w]
```

where `A` and `r_r` are reservoir contributions, `D` and `r_w` are the well
system, and `B`/`C` couple reservoir cells and well unknowns. Multisegment
wells use the same interface with a larger sparse well block.

The coupling crosses the reservoir boundary in three steps:

1. During TPFA assembly, `addReservoirSourceTerms()` inserts each connection
   rate and its automatic-differentiation derivatives into the reservoir cell
   residual and diagonal block.
2. Before the linear solve, `WellConnectionAuxiliaryModule::linearize()` adds
   the eliminated matrix contribution when configured and always applies the
   eliminated residual contribution. Conceptually this produces
   `A - C^T D^-1 B` and `r_r - C^T D^-1 r_w`.
3. After the reservoir solve, `postSolve()` gathers `delta_r` at each well's
   perforated cells and recovers
   `delta_w = D^-1 (r_w - B delta_r)`, then updates `WellState`.

This protocol is mathematically cohesive, but its interface is distributed
among TPFA source assembly, an auxiliary-module hook, and a post-solve hook.

### Convergence

`getWellConvergence()` gathers per-well residual failures across ranks. When
requested by the global nonlinear solver it also adds two controller-level
conditions:

- the latest `well_group_control_changed` flag is reported as violated
  well/group targets;
- a network that must continue balancing forces another Newton iteration.

The result deliberately combines equation convergence with discrete control
stability and network fixed-point progress. This is useful operationally but
makes “well convergence” mean more than convergence of the assembled well
equations.

### Successful timestep and rollback

At a successful timestep, `timeStepSucceeded()` updates well-throughput and
connection histories, injection multipliers, potentials, well-test/economic
state, productivity indices, and selected group modes. It then calls
`commitWGState()`, which copies the active WG state to the accepted copy and
commits network node and branch data.

If a later timestep attempt fails, the reservoir state is rolled back by the
reservoir lifecycle and the next well `beginTimeStep()` restores the accepted
WG/network copies. No explicit `BlackoilWellModel::updateFailed()` mutation is
needed for that copy; retry entry performs the reset.

The transaction boundary is nevertheless incomplete as an abstraction:

- `WGState` is copied as one unit, but network state is committed separately.
- The NUPCOL snapshot is a third copy with a policy-specific lifetime.
- Well computational objects contain mutable state that is not represented by
  `WGState` and is re-established through creation and initialization calls.
- Guide-rate data, filter-cake/injection histories, open/close bookkeeping,
  and caches have their own lifecycle rules.
- `guideRate_`, the auto-choke cache, and `last_glift_opt_time_` are serialized
  but are not members of the accepted WG snapshot. Their retry behavior comes
  from recomputation or algorithm-specific rules rather than a uniform restore.
- `BlackoilWellModelNetwork::doPreStepRebalance()` explicitly contains a TODO
  asking whether group and network state should be backed up when its local
  solve does not converge.

## Well subsystem

### Inputs from the reservoir

Well code reads reservoir data primarily through `Simulator`:

- intensive quantities and fluid state at perforated cells;
- mobilities, PVT quantities, component mappings, densities, and
  transmissibility multipliers;
- current timestep size, simulation time, and Newton iteration context;
- average formation-volume factors used for convergence scaling;
- the current reservoir Newton update after the linear solve.

The set is larger than the signature of most `BlackoilWellModel` entry points
suggests because a well receives the full `Simulator` or can reach it through
its owner.

### Persistent state versus equation state

`SingleWellState` contains the serializable results and controller state,
including:

- status, BHP, THP, temperature, efficiency, and control modes;
- surface, reservoir, previous, connection, mixing, energy, and fracture
  rates;
- perforation and segment state;
- potentials, productivity indices, and implicit IPR coefficients;
- primary-variable values;
- group target and fallback target;
- retained network THP limit;
- ALQ value and gas-lift increment/decrement counters;
- event and WELDRAW state.

The `WellInterface` object separately owns or tracks:

- the concrete primary-variable representation used by the equations;
- B/C/D matrices and the well residual;
- connection-rate automatic-differentiation values;
- operability and solvability state;
- dynamic THP limit and control/status switching machinery;
- local convergence, reopening, and switching counters/logs.

Several methods synchronize the two representations, for example
`updatePrimaryVariables()`, `updateWellState()`, and network initialization.
This two-way synchronization is one source of non-obvious state flow.

### Outputs to the reservoir and solver

The well subsystem exposes four numerically distinct products:

1. connection source terms for reservoir residual/Jacobian assembly;
2. well convergence failures and controller-stability flags;
3. eliminated well matrix/residual contributions for the reservoir linear
   system and optional preconditioner/GPU representations;
4. updated well unknowns reconstructed after the reservoir solve.

These are a coherent “well coupling result” conceptually, but no object in the
current interface groups them or states their validity lifetime. They remain
valid only relative to the reservoir state and controls at which the last well
assembly was performed.

## Group-control subsystem

### Aggregation and target allocation

`updateAndCommunicateGroupData()` is the main bridge between per-well state
and group control. It:

1. establishes which wells are globally under group control;
2. chooses or refreshes the NUPCOL snapshot;
3. updates counts of group-controlled wells;
4. recursively calculates injection, production, REIN, VREP, reduction, and
   network-leaf rates;
5. updates active well rates using the selected current/frozen state;
6. performs global communication of well and group rates;
7. optionally derives a `group_target` and producer fallback target for each
   GRUP-controlled well.

Target allocation depends on the schedule group hierarchy, active group
control modes, guide rates, efficiency factors, rate-conversion coefficients,
and reduction rates. If no usable group target can be obtained, a GRUP well
can be changed to BHP control.

NUPCOL is not simply a convergence limit. It defines a data snapshot used to
freeze selected well/group quantities after the early Newton iterations.
Some quantities, such as group target reductions, intentionally continue to
be updated. There is also a relative-rate-change exception that can refresh
the NUPCOL snapshot for REIN/VREP behavior. Consequently “current group data”
can mix values derived from current and frozen well state.

### Constraint switching

`updateWellControls()` recursively checks higher-level and individual group
constraints and then well constraints. When a group or well changes control,
it calls `updateAndCommunicate()` so all ranks see the new modes and targets.
That helper performs group aggregation, updates GRUP wells from their targets,
and aggregates again.

The control algorithm therefore has nested feedback even before the network
loop:

```text
check group constraint
  -> mutate GroupState control mode
  -> globally recompute rates/targets
  -> update GRUP-controlled WellState
  -> globally recompute rates/targets again
  -> check the next constraint
```

Control-change return values are reduced across MPI ranks, but the detailed
result remains encoded in mutated state and switching logs.

### `GroupStateHelper` as an implicit interface

`GroupStateHelper` provides schedule, summary, group hierarchy, targets,
rates, guide rates, conversion helpers, logging, and optional reservoir-
coupling services. It holds pointers/references to the current `WellState` and
`GroupState`.

RAII guards such as `pushWellState()` and `pushGroupState()` temporarily
retarget the helper to alternate objects—for example the NUPCOL state—then
restore the previous target. This avoids threading many state parameters
through recursive functions, but it makes the actual state read by a helper
call depend on ambient guard scope. A method signature taking only
`GroupStateHelper&` does not reveal whether it operates on active, NUPCOL, or
temporary copied state.

## Production-network subsystem

### Inputs and outputs

The extended network consumes:

- schedule network topology, terminal pressures, branch properties, and
  NETBALAN settings;
- globally aggregated production/injection rates at network leaf groups;
- VFP production properties and unit system;
- cached auto-choke group THP values;
- current well and group state.

It produces:

- node pressures and branch flow data;
- a maximum pressure-change imbalance for convergence;
- dynamic THP limits on prediction-mode producer wells with VFP tables;
- auto-choke group THP stored in `GroupState`;
- output node/branch records.

### Network fixed point

`BlackoilWellModelNetwork::update()` first solves auto-choke group THP when
needed. It then runs network sub-iterations:

```text
aggregate leaf rates
  -> compute undamped node pressures and branch flows
  -> compare with previous node pressures
  -> damp and cap each pressure update
  -> impose leaf pressure as dynamic well THP limit
  -> if not converged:
       locally prepare/solve affected producer wells
       aggregate group and leaf rates again
       repeat
```

The globally reduced maximum node-pressure change is compared with the
NETBALAN pressure tolerance. After sufficiently many outer iterations the
tolerance can be relaxed by a factor of ten. Balancing is scheduled either at
timestep start or within the NUPCOL window; unsupported time-interval mode is
not implemented here.

The network is thus coupled by fixed-point iteration rather than by adding
node-pressure unknowns and network equations to the reservoir/well Jacobian.
That choice keeps the reservoir linear system smaller but distributes
convergence across global Newton, network outer iterations, network
sub-iterations, and local well iterations.

### Auto-choke groups

For an auto-choke node, `computeWellGroupThp()` derives the applicable group
target and root-finds a common THP. Each mismatch evaluation applies a trial
dynamic THP and solves the group's well equations to calculate group rate.
The accepted THP is cached, stored in `GroupState`, and applied to eligible
wells.

This is an especially deep state interaction: a network calculation performs
well solves, mutates dynamic constraints and well state during root finding,
and feeds a group-level result back into both group and well representations.

### Current, accepted, and reported pressure

The network owns tentative and accepted maps for node pressures and branch
data. They participate in `commitWGState()`/`resetWGState()` through separate
network calls.

For output, `assignNodeAndBranchValues()` exports the stored, damped runtime
values and also calls the pressure computation once more to produce a
“converged” undamped pressure/branch result. Output is therefore not a pure
serialization of the runtime maps; it performs a derived network calculation.

## Gas-lift subsystem

Gas lift is invoked after network pressure update and before the final local
well preparation in each network outer iteration.

`maybeDoGasLiftOptimize()` first checks schedule activation and the minimum
wait-time policy. When network balancing may have changed well THP, it
recomputes well potentials using the current network node pressures. It then
constructs temporary `GasLiftGroupInfo` and per-well optimization objects.

The two stages are:

1. **Stage 1: per-well optimization.** Eligible local wells add or remove ALQ
   while observing group limits. Ranks execute this stage serially and
   broadcast changed group oil, gas, water, and ALQ totals because later wells
   need the latest allocation.
2. **Stage 2: group reallocation.** The group hierarchy redistributes or
   removes increments using production gradients and constraints. This updates
   `SingleWellState::alq_state` and well potentials.

Most optimization structures are temporary. Persistent gas-lift state is
split between:

- `BlackoilWellModelGasLift::last_glift_opt_time_`;
- each well's serialized `ALQState` (current/default ALQ and oscillation
  counters);
- well potentials in `SingleWellState`;
- the general guide-rate state owned by the well model.

If any rank changes a well, the change is reduced globally. The caller treats
that as a reason to update guide rates and perform another network outer
iteration when network balancing is active. This makes the practical data
flow:

```text
network pressure -> potential -> ALQ -> well rate -> group leaf rate
       ^                                             |
       +---------------------------------------------+
```

The optimizer does not return a structured allocation delta. Its result is a
boolean plus mutation of well potentials, ALQ state, and temporary/group
accounting structures.

## Guide rates: the coupling layer between group control and gas lift

Guide rates deserve explicit mention even though they are not one of the four
named subsystems. They influence how group targets are distributed to wells
and groups. They are updated at timestep start, when their schedule policy
requires it, and after ALQ changes. Gas lift changes potentials; potentials
can change guide rates; guide rates change group target allocation; target
allocation changes well rates; and those rates drive the network and gas-lift
constraints.

`GuideRate` is serialized as part of `BlackoilWellModelGeneric`, while
`GuideRateHandler` coordinates updates using schedule, summary, group, well,
and optional reservoir-coupling data. This is another case where the durable
state and its algorithm/service object are separate.

## Parallel data flow

Parallel behavior cannot be treated as an implementation detail because it
shapes the state types and call ordering.

| Data or decision | Parallel behavior |
| --- | --- |
| Distributed-well equation blocks | Perforation contributions to well residual and `D` are summed over the well communicator; `B`/`C` use distributed well structures |
| Per-well mixing and rates | Selected quantities are summed over the parallel well communicator |
| Global well/group rates | `WellState::communicateGroupRates()` and `GroupState::communicate_rates()` aggregate maps/vectors across the grid communicator |
| Control changes | Boolean/count changes are reduced so all ranks repeat compatible control logic |
| Well convergence | Local reports are gathered into a global `ConvergenceReport` |
| Network convergence | Maximum node-pressure imbalance is reduced globally |
| Network activity | A global maximum decides whether any rank activates the network |
| Gas-lift stage 1 | Ranks take turns optimizing local wells; updated group totals are broadcast after each rank |
| Gas-lift changed flag | Number of changed wells is globally summed |

`WellState` reflects two views at once. Its `wells_` container holds wells
present on the rank, and a distributed well may appear on multiple ranks.
Other maps such as global well information/rates have entries for all wells.
Code must therefore distinguish local well index, global well name, well
communicator, and grid communicator.

The multiple synchronization scopes are a major reason that apparently local
functions have strict ordering constraints. Replacing shared mutation with
return values alone will not simplify the design unless those return values
also make ownership and collective scope explicit.

## Restart, output, and serialization

The well model serializes considerably more than `WGState`:

- report-step flags and well/perforation bookkeeping;
- guide-rate state;
- current and accepted network pressure/branch state and auto-choke cache;
- injection-multiplier and filter-cake state;
- active, accepted, and NUPCOL `WGState` copies;
- switched-group and offending-well sets;
- gas-lift optimization time.

Within `SingleWellState`, network THP, ALQ, group targets, primary variables,
perforations, and segments are serialized. `GroupState` serializes its main
rate/control/reduction maps, network leaf rates, group THP, constraint targets,
and well counts. Some construction and derived/global lookup information in
`WellState` is deliberately reconstructed instead of serialized.

Output has separate extraction paths:

- `wellData()` reports `WellState`, then enriches it with guide rates,
  targets, dynamic status, tracers/species, and shut-connection data;
- `groupAndNetworkData()` maps `GroupState` controls and guide rates into
  report data and asks the network for node/branch values;
- the writer collects rank-local objects onto the I/O rank.

This means the report schema is not the internal state schema. Output assembly
is another integration layer with schedule and derived calculations.

## Why the data flow is hard to understand

### 1. Mutation is the primary API

Calls such as `updateAndCommunicateGroupData()`, `updateWellControls()`,
`network.update()`, `prepareWellBeforeAssembling()`, and gas-lift optimization
all mutate overlapping state. A returned `bool` says that something changed,
but the changed fields and their validity are implicit.

### 2. Temporal meaning is encoded in member names and call order

Active, accepted, previous, NUPCOL, initial, current, and report values exist,
but they are not represented by a common timestep-state abstraction. The
caller must know, for example, whether “previous well state” means the accepted
timestep copy or the NUPCOL-frozen copy.

### 3. Derived and authoritative data share the same objects

`WellState` mixes primary/controller state with potentials, IPR values,
targets, cached previous rates, and reported quantities. `GroupState` mixes
control decisions with aggregates that can be recomputed. Network state mixes
iterated pressure with output branch data. This obscures which values must be
rolled back, serialized, invalidated, or recalculated.

### 4. Service-location hides dependencies

The full `Simulator` and a broad `GroupStateHelper` give algorithms access to
reservoir properties, schedule, summary state, communicators, logging, guide
rates, and mutable state. Signatures therefore understate both read and write
sets.

### 5. Several nonlinear loops overlap

A single global Newton iteration can contain group-constraint loops, network
outer loops, network pressure sub-iterations, auto-choke root finding, local
well nonlinear iterations, control/status switching, and two-stage gas-lift
optimization. “Iteration” without a qualifier is ambiguous.

### 6. Numerical coupling and operational policy are interleaved

Well B/C/D assembly sits beside economic closure, NUPCOL freezing, guide-rate
policy, control switching, network balancing schedule, and gas-lift wait-time
policy. These concerns legitimately interact, but their state transitions are
not separated at the interface.

### 7. The facade has too many consumers

Reservoir assembly, nonlinear convergence, linear solvers/preconditioners,
output, restart, actions, well testing, and reservoir coupling all call into
`BlackoilWellModel`. Accommodating these consumers has produced multiple
specialized accessors and alternative extraction paths.

## Candidate interface simplifications

The following sequence is intended to improve observability before attempting
large algorithmic changes.

### 1. Define explicit state categories

Document and eventually type the state as four categories:

```text
AcceptedFacilityState
    physical/controller state committed at timestep success

FacilityIterateState
    tentative well, group, network, and ALQ state for the current attempt

FacilityPolicySnapshot
    NUPCOL-frozen rates/targets and other iteration-policy values

FacilityDerivedData
    aggregates, potentials, guide rates, cell rates, equation blocks, output caches
```

The immediate gain would be a reviewable rule for serialization, rollback,
and invalidation of every field.

### 2. Make one assembly epoch explicit

Introduce a result or token representing a well assembly at a particular
reservoir/control epoch. It could own or reference:

- connection sources;
- well B/C/D/residual blocks;
- convergence result;
- cell-to-well mapping;
- the state/version against which recovery is valid.

Reservoir source insertion, Schur elimination, and post-solve recovery would
then be visibly related operations rather than independent callbacks.

### 3. Return structured controller deltas

Replace boolean-only results incrementally with records such as:

```text
ControlUpdateResult
    changed wells/groups, old/new modes, targets requiring recomputation

NetworkUpdateResult
    node pressures, branch data, imposed THP limits, imbalance, continue reason

GasLiftUpdateResult
    ALQ deltas, potential deltas, affected groups, guide-rate invalidation
```

The implementation may still apply these deltas initially, but exposing them
would make dependencies and convergence reasons testable.

### 4. Split group state from group aggregates

Keep durable control modes, GPMAINT/controller history, and auto-choke state in
one object. Put recomputable rates, reductions, counts, and network-leaf totals
in an explicitly versioned `GroupAggregates` object. Target allocation should
state whether it consumes current or NUPCOL aggregates.

### 5. Narrow `GroupStateHelper`

Separate at least:

- a read-only schedule/group hierarchy view;
- a rate/target calculation service with explicit input state;
- a communication service;
- a controller mutation service;
- logging and reservoir-coupling adapters.

As a first low-risk step, functions using alternate NUPCOL state could accept
that state explicitly instead of relying on `pushWellState()`/
`pushGroupState()` pointer retargeting.

### 6. Give the network a side-effect-light solve boundary

A network solve could consume topology plus immutable leaf rates and return
candidate node pressures, branch flows, and imbalance. Applying damping,
committing pressure state, imposing well THP limits, and resolving wells could
be separate named steps. Auto-choke root finding should ideally evaluate trial
well solutions in scratch state rather than mutating the active state through
each mismatch evaluation.

### 7. Isolate gas-lift allocation from state application

Gas-lift stages could return an allocation plan based on explicit well/group
snapshots. Applying ALQ, updating potentials, communicating group totals, and
invalidating guide rates would then be visible orchestration steps. This would
also make the currently serialized versus temporary gas-lift state boundary
easier to test.

### 8. Name every loop level

Use distinct result types and counters for:

- global reservoir Newton iteration;
- facility/network outer fixed point;
- network pressure sub-iteration;
- group-constraint iteration;
- local well nonlinear iteration;
- auto-choke root iteration;
- gas-lift stage iteration.

This is partly terminology, but it would remove ambiguity from diagnostics,
parameters, convergence reports, and future interfaces.

## Suggested implementation order for follow-up work

1. Build a state-field inventory marking each member as authoritative or
   derived, its owner, commit/rollback rule, restart behavior, and MPI scope.
2. Add focused trace tests around one group-controlled producer, one active
   network leaf, and one gas-lift well. Record control modes, target version,
   node pressure, THP limit, ALQ, and rates at every loop boundary.
3. Add debug-only generation counters to well assembly, group aggregates,
   guide rates, and network leaf rates. Assert that consumers use a compatible
   generation.
4. Extract a read-only `ReservoirWellView` from `Simulator` for perforation
   properties and iteration metadata.
5. Make group aggregation return an explicit `GroupAggregates` value while
   retaining the current mutation path as an adapter.
6. Make network and gas-lift updates return structured result records.
7. Only after those boundaries are observable, reconsider whether network
   pressure or facility controls should remain fixed-point coupled or move
   closer to the global nonlinear system.

## Source map

The main files inspected were:

- `opm/simulators/wells/BlackoilWellModel.hpp`
- `opm/simulators/wells/BlackoilWellModel_impl.hpp`
- `opm/simulators/wells/BlackoilWellModelGeneric.hpp`
- `opm/simulators/wells/BlackoilWellModelGeneric.cpp`
- `opm/simulators/wells/WellInterface.hpp`
- `opm/simulators/wells/WellInterface_impl.hpp`
- `opm/simulators/wells/StandardWell_impl.hpp`
- `opm/simulators/wells/StandardWellEquations.hpp`
- `opm/simulators/wells/StandardWellEquations.cpp`
- `opm/simulators/wells/MultisegmentWell_impl.hpp`
- `opm/simulators/wells/MultisegmentWellEquations.cpp`
- `opm/simulators/wells/WellConnectionAuxiliaryModule.hpp`
- `opm/simulators/wells/WellState.hpp`
- `opm/simulators/wells/SingleWellState.hpp`
- `opm/simulators/wells/WGState.hpp`
- `opm/simulators/wells/GroupState.hpp`
- `opm/simulators/wells/GroupStateHelper.hpp`
- `opm/simulators/wells/GroupStateHelper.cpp`
- `opm/simulators/wells/BlackoilWellModelGuideRates.cpp`
- `opm/simulators/wells/GuideRateHandler.hpp`
- `opm/simulators/wells/GuideRateHandler.cpp`
- `opm/simulators/wells/BlackoilWellModelNetwork.hpp`
- `opm/simulators/wells/BlackoilWellModelNetwork_impl.hpp`
- `opm/simulators/wells/BlackoilWellModelNetworkGeneric.hpp`
- `opm/simulators/wells/BlackoilWellModelNetworkGeneric.cpp`
- `opm/simulators/wells/BlackoilWellModelNetworkPressureComputation.hpp`
- `opm/simulators/wells/BlackoilWellModelGasLift.hpp`
- `opm/simulators/wells/BlackoilWellModelGasLift_impl.hpp`
- `opm/simulators/wells/GasLiftGroupInfo.hpp`
- `opm/simulators/wells/GasLiftGroupInfo.cpp`
- `opm/simulators/wells/GasLiftSingleWell_impl.hpp`
- `opm/simulators/wells/GasLiftStage2.hpp`
- `opm/simulators/wells/GasLiftStage2.cpp`
- `opm/simulators/wells/ALQState.hpp`
- `opm/simulators/flow/FlowProblem.hpp`
- `opm/simulators/flow/NonlinearSystemBlackOilReservoir_impl.hpp`
- `opm/simulators/flow/EclWriter.hpp`
