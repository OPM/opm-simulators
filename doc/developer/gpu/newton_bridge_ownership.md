# Gas–water–energy Newton bridge ownership

`GpuFlowGasWaterEnergyBridge<CpuTypeTag, DeviceTypeTag>` is the only owner of
resident solution state. The dispatcher owns the bridge, and assembly borrows
views. No view extends an allocation's lifetime.

The bridge retains current and previous typed primary-variable buffers,
current and previous IQ buffers, a transactional primary-variable scratch
buffer, current and scratch byte switch histories, current and unrelaxed
previous correction vectors, and four update-status words. It also retains the
fluid system, material/thermal laws, rock data and cell volumes. Typed PV
conversion preserves every meaning flag, PVT region and capillary metadata;
CPU and GPU object layouts are never assumed equal.

## Generations and publication

Each time slot has a primary-variable generation, a host-mirror generation and
an IQ generation. Zero denotes invalid data. Importing explicitly edited CPU
state increments the resident generation, makes its host mirror current, and
invalidates IQs. Device Newton writes only scratch state. After a successful
compact-status check, commit swaps scratch into current, increments the PV
generation, invalidates IQs and retains the previous accepted time slot. Failed
kernels never publish their partial scratch state.

Property completion associates IQs with the resident PV generation. Host
materialization only updates the mirror generation: it neither invalidates
resident state nor requests a later upload. Callers must explicitly notify the
bridge after CPU mutations; obtaining a mutable reference is not a mutation.
The materialization destination is the canonical model solution for that slot,
not an unrelated temporary vector.

On accepted timesteps, device copies advance current PVs and IQs into history.
On rejected attempts, previous PVs are copied into current, given a new
generation, and current IQs must be reevaluated. Correction history is reset
where CPU `dx_old_` resets. Switch history deliberately survives both transitions,
matching the CPU Newton method's `wasSwitched_` lifetime.

## Streams and lifetime

The existing solver and assembly execute on the default stream. Properties,
PV updates, and lifecycle copies use the bridge stream. Before a solve, an
event orders default-stream work after bridge writes. After the solve, a second
event orders the bridge update after solver/assembly completion. Property-ready
events order assembly consumers. No device-wide synchronization is used for
these transitions. The compact status download waits for the bridge stream
before checking failure or publishing scratch state.

Replacing or destroying owners captures completion of default-stream readers
and waits for bridge writers before freeing allocations. Host import waits at
its explicit compatibility boundary so the retained staging vector cannot be
rewritten while its previous upload is in flight.

## Explicit transfer boundaries

* Initialization/reinitialization uploads static owners, history seeds and the
  initial primary-variable states, and allocates bridge buffers.
* Explicit host-state import is permitted after initialization, restart loading
  or a CPU fallback mutation. Ordinary resident Newton updates do not import.
* Primary-variable and IQ materialization serve CPU convergence, diagnostics,
  output and compatibility consumers. Repeated PV reads of a current mirror do
  not transfer again.
* Correction materialization is only for explicitly requested diagnostics or
  shadow validation. Ordinary resident solves do not download corrections.
* Switch-history imports/downloads serve CPU fallback or shadow validation.
  Activation also imports CPU `dx_old_` into the existing previous-correction
  buffer so SOR history survives a CPU fallback. This exceptional upload is
  counted separately from ordinary resident updates.
* Every update downloads four 32-bit status words: failing cell plus one,
  failure reason, switched-cell count and a reserved word.

Bridge transfer counters distinguish PV uploads/downloads, correction downloads,
IQ downloads, successful updates and owned-buffer allocations. Scoped low-level
accounting records every allocation and upload through `GpuBuffer`, `GpuVector`
and GPU smart pointers during bridge initialization. Nested static-data scopes
measure the fluid system, problem data and volumes separately from evolving PVs
and IQs. The hooks are inactive outside their constructing host thread and do
not log individual operations. Static upload batches additionally count owner
initializations. Assembly's independent parameter-buffer construction remains
outside bridge allocation counters. Host solver compatibility calls count their
temporary correction allocation separately from persistent bridge allocations.
