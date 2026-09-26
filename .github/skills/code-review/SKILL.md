---
name: code-review
description: >-
  The domain checklist for reviewing pull requests in OPM/opm-simulators,
  OPM/opm-common, OPM/opm-grid and OPM/opm-upscaling: severity/confidence
  rubric, blocking-issue criteria, path-to-section routing, domain passes
  (restart/serialization, deck/parser/keywords, wells/groups/networks,
  numerics, determinism, parallel/MPI, gpuistl), cross-cutting passes (C++
  hygiene, headers/build/CMake, performance, tests, style), and CI/merge
  conventions. Apply this whenever reviewing a pull request against these
  repositories: what to look for, and how to grade what you find.
---

# OPM Pull Request Review Checklist

This encodes the review standards applied in the OPM repositories: what counts as a defect, how to grade it, and which passes apply to which parts of the codebase. It says nothing about contributor-facing build/format/test mechanics, which live in `CONTRIBUTING.md` in `opm-simulators`.

Ownership in this project is broad and shared — most areas have several people who could weigh in, not one designated owner. Don't name an individual to route a finding to; when something is beyond what you can verify, say so as a Question for the PR thread.

The checks below name real symbols and files. Prefer verifying a check with a `grep`/`Read` over reasoning about it from the diff alone — a verified finding is worth five speculative ones.

## Severity and confidence

Every finding gets both. They are independent axes: a high-confidence nit and a low-confidence blocker are both possible, and the reader triages on the pair.

**Severity**

- **Blocker** — merge must not happen. Restricted to the classes below, plus undefined behaviour, data races, lifetime bugs, and non-determinism in result-critical paths.
- **Should-fix** — real defect risk in edge cases, weak error handling on invalid decks, likely hot-path regression without evidence, API change without a migration note.
- **Nit** — clarity, naming, hygiene. Goes in the Nits section; don't repeat the label inline.
- **Question** — you need author knowledge to grade it. Anything you cannot ground in the diff goes here, not in Blocking.

Do **not** promote a finding to Blocker merely because the PR text does not state intent. Undocumented intent is a Question. Reserve `CHANGES_REQUESTED` for architectural objections, file-format breakage, or missing tests; everything else is inline comments plus `COMMENTED`, or `APPROVED` with nits.

**Should-fix or Question?** These two collide constantly, because the commonest thing you want to say is "this looks wrong — was it intended?". Decide on one test: *can you describe a concrete path to failure in the code as written?*

- You can name the input, deck feature, phase configuration, rank count or state that makes this behave wrongly → **Should-fix**, whatever the PR text says about intent. Put the intent question in that entry's own text ("this shuts wells that previously stayed open under X — intended?"). Do not also list it under Questions; one finding, one section.
- Grading it needs knowledge you do not have — what testing indicated a value, whether a deck can reach the branch, what the author meant → **Question**. Nothing under Questions gets a severity, because you are not asserting a defect.
- The only thing wrong is that the PR does not explain itself → **Question**, per the rule above.

If the answer changes once the author replies, that is expected and fine: a Question can become a Should-fix later. Guessing in order to avoid asking is what you must not do.

**Confidence**

- **High** — you read the code path, or verified it with a grep/`Read` against the checkout. Say what you checked.
- **Medium** — strongly implied by the diff, not independently verified.
- **Low** — pattern match only. Low-confidence items are phrased as Questions and never appear under Blocking.

## Routing by touched path

Apply the listed sections hard; skim the rest. Match review depth to the change — a one-line fix does not get the whole checklist.

| Path | Sections |
|---|---|
| `opm/output/eclipse/`, `opm/io/eclipse/`, `*RestartIO*`, `*Summary*` | Restart, Tests |
| `opm/input/eclipse/`, `opm/input/eclipse/share/keywords/` | Deck, Restart, Tests |
| `opm/simulators/wells/`, `*Group*`, `*Network*` | Wells, Determinism, Parallel, Numerics |
| `opm/simulators/linalg/gpuistl/`, `*.cu` | GPU |
| `opm/simulators/linalg/` (CPU), `*Preconditioner*`, `*ISTLSolver*` | Numerics, Determinism, Performance |
| `opm/simulators/timestepping/`, `opm/models/`, `opm/material/` | Numerics |
| `opm/grid/cpgrid/`, `opm/grid/common/` (partitioning, load balance) | Determinism, Parallel, Performance |
| `opm/upscaling/`, `opm/porsol/`, `opm/elasticity/` | Numerics, Build |
| `CMakeLists*.cmake`, `CMakeLists.txt`, `jenkins/` | Build |
| Anything with new/removed files | Build (file lists), Style |

## Blocking issues

These are the objections that actually stop merges. Raise them before any line-level comment.

- **Restart / summary / INIT file layout changed.** Any change to a `VectorItems` enumerator (`VI::intehead::`, `VI::XGroup::`, `VI::SWell::`), array ordering, or an integer written to the restart file is a compatibility change and must be called out explicitly with a reference-solution plan.
- **PR does too many things.** Refactor + semantics + formatting in one branch. Ask for a stacked review, or for a split into ordered, individually reviewable pieces, and for semantics-changing hunks (e.g. `epsilon` → `0.0`) to be isolated.
- **New non-trivial logic with no unit test.** Recommend a unit test covering the new logic and name the file it belongs in (`tests/*Tests.cpp`). Non-trivial logic is not merged untested.
- **Bad layering / coupling.** A type from a heavier layer taken as a constructor parameter purely for convenience. Propose a callback or a narrower interface — the pattern to follow is `SummaryConfig`, which takes `std::function<GridDims(const std::string&)>` rather than an `EclipseGrid` (`SummaryConfig.hpp`).
- **Behaviour change presented as a refactor.** A pure refactor must produce equal output up to given precision. If the PR is described as refactor-only but a hunk changes a comparison, an ordering, a default or a tolerance, it is a behaviour change and has to be labelled as one.
- **Non-determinism introduced in a result-critical path.** See the Determinism pass.

## Domain passes

### Restart, serialization, backward compatibility
- Serialize only dynamic state that is *not* recreated at report-step start and *not* derivable from the primary solution.
- Every new serializable type needs `serializeOp()`, a `serializationTestObject()` built from **non-default** values, and an entry in `test_Serialization` / `test_RestartSerialization`.
- Restart items used to re-initialise UDAs must be stored in **output** units, or `eval_uda()` converts twice.
- New restart-relevant vectors must be added to `requiredRestartVectors()` in `Summary.cpp`.
- Don't reconstruct state from a restart file that never carried it; guard reconstruction on whether the deck defines the feature at all.
- Input and output file conventions must stay consistent (`IOConfig::consistentFileFlags()`); derive unified-vs-`.X000N` from how the file was written, never assume.

### Deck, parser, keywords (opm-common)
- Use compiled items `record.getItem<ParserKeywords::KW::ITEM>()`, never string lookups; check `defaultApplied(0)` / `hasValue(0)`; use `getSIDouble(0)` / `getTrimmedString(0)`.
- Validate before storing: get a valid parsed value first, then assign to the member.
- Deck errors are `OpmInputError`, or `OpmLog::warning(OpmInputError::format(msg, location))` — never hand-formatted file/line strings. `assert` only for what the parser already guarantees.
- Changing what a *defaulted* item resolves to is a design change; require confirmation against the product owner.
- Name matching uses `shmatch()` (`opm/common/utility/shmatch.hpp`) or the existing `wellNames(pattern)` helper.
- New time-dependent `ScheduleState` members use `map_member<K,V>`, not a raw `unordered_map`; don't store a whole `UnitSystem` per record.
- Newly supported items must move out of `PartiallySupportedFlowKeywords.cpp`.
- No speculative API: don't add a `ScheduleEvents` bit, member or public function that nothing in the PR uses.

### Wells, groups, networks
- A fix in `StandardWell_impl.hpp` almost always needs the mirror in `MultisegmentWell_impl.hpp`. Ask every time.
- Control-switch predicates must be justified in both directions (injector *and* producer), stating which control is more restrictive.
- Flags must be true at their initial value: name and default a new boolean so the safe/conservative state is the default, not the optimistic one.
- Defaults must not fall into a fallback path by accident. Check each new default against every threshold it is compared with.
- Network tolerances come from `Network::Balance` (a percentage of the group rate target), not a hard-coded `1e-4`; respect `NETBALAN` and `nupcol`.
- Static deck config lives in `Schedule`; dynamic state in `WellState` / `GroupState` / `WellTestState`. Don't mix.
- Never index well state by `Well::seqIndex()` — it only matches in serial. Use `wellState().index(name)`.
- Connections are grid-distributed; segments are not. Don't assume they distribute together.
- Anything that makes a well more likely to end up SHUT is a behaviour change. The project's default position is that a well stays OPEN as long as it can operate, so a change that shuts wells earlier needs explicit justification.
- Guard derivative index blocks (`WFrac`, `GFrac`) on `FluidSystem::phaseIsActive()`; `canonicalToActivePhaseIdx()` can throw and must not be called inside an OpenMP region.

### Numerics and physics
- Question every changed tolerance, epsilon, safety factor, iteration limit and default, and ask what testing indicated the new value. A parameter change alters the equation being solved and belongs in its own PR.
- No magic dimensional constants. Use a named constant in `Opm::unit::` terms (`5.0 * unit::barsa`, `unit::day`). Duplicated epsilons for the same physical situation must be hoisted to one header.
- Check that a new tolerance can actually fire by evaluating it numerically — multiplying two already-small quantities together can drop the threshold far enough that it never triggers.
- Assigning a plain constant into an `Evaluation` destroys derivatives.
- Energy/mass coupling must use the same density, the same upwind state and the same derivative source as the mass balance.
- Iterate active phases (`activeToCanonicalPhaseIdx`, `canonicalToActiveCompIdx`), never all three canonical phases; phase-indexed access assuming three phases is out-of-bounds in a gas/water run.
- Guard near-zero denominators and floating-point equality; fixed epsilons degenerate to zero under `float` instantiation.
- Prefer `std::optional` to sentinel values; distinguish "legitimately absent" (return sentinel) from "structural error" (`OPM_THROW`).
- Watch unsigned subtraction; `integer / integer` is integer division.
- Template on `Scalar`; guard `double`-only paths with `if constexpr (std::is_same_v<Scalar, double>)`.

### Determinism and reproducibility
Results must not depend on container hashing or allocation addresses.

- Flag iteration over `unordered_map` / `unordered_set` / pointer-keyed containers where the loop body accumulates into a result, orders output, or decides control. Iterate a sorted key vector instead.
- `std::stable_sort` where insertion order matters; never sort by pointer address or by a hash.
- Changing the order or grouping of a floating-point accumulation changes the answer. Reductions over cells, connections or ranks need a fixed order, and a new OpenMP `reduction` clause on a previously serial loop is a behaviour change — say so.
- Group/network decisions must be rank-invariant: all ranks see the same information and reach the same decision.
- Results must not depend on `OMP_NUM_THREADS` or on the domain-decomposition partition unless the PR states that they do and why.
- A "no change" claim needs runtime evidence, not an argument. Ask for the deck name and iteration counts (do not re-run the case yourself).

### Parallel and MPI
- No new `MPI_COMM_WORLD` or default-constructed collective communication. No `MPI_Abort` — return `EXIT_FAILURE` and let `Main::~Main()` finalize.
- Exceptions must propagate to a collective handler: check `OPM_BEGIN_PARALLEL_TRY_CATCH` / `OPM_END_PARALLEL_TRY_CATCH` placement, and `OPM_DEFLOG_THROW` with `DeferredLogger` in well code. A bare throw on some ranks deadlocks.
- Collectives inside per-well loops whose iteration order may differ between ranks are a deadlock risk.
- Never communicate non-POD types (`std::pair<std::string,double>`) through `allgather`.
- Mask on `Dune::InteriorEntity` in every accumulation; use `collectToIORank_.isIORank()`, not `comm.rank() == 0`.

### GPU / gpuistl
For the class-porting pattern (`GpuBuffer`/`GpuView` ownership, `OPM_HOST_DEVICE` decoration, `copy_to_gpu`/`make_view`, testing against a CPU reference) see `opm-simulators/doc/developer/gpu/porting_to_gpu.md` — verify against its checklist rather than re-deriving the rules here. What that doc does not cover, and still needs a reviewer's eye:
- Ownership beyond `GpuBuffer`: every stream, event, graph, and cuSPARSE/cuBLAS handle or descriptor is RAII-owned. A throwing `*_SAFE_CALL` must not leak.
- Wrap every GPU API call in the matching `OPM_GPU_SAFE_CALL` / `OPM_CUSPARSE_SAFE_CALL` / `OPM_CUBLAS_SAFE_CALL` / `OPM_HYPRE_SAFE_CALL`.
- A class creating its own stream must destroy *and* synchronise it. Async functions get an `Async` suffix; the default stream is a named constant, never a literal `0`.
- Portability: guard cuSPARSE generic-API use by CUDA version; every capability `#if` needs an `#else` that errors or falls back, never silent degradation; add `<type_traits>` wherever `is_same_v` appears.
New GPU types go into `is_gpu_type` (`gpuistl/detail/gpu_type_detection.hpp`) / `is_gpu_operator_v` (`linalg/is_gpu_operator.hpp`); CPU-only classes must still compile with GPU support on; GPU solvers must match `ISTLSolver` semantics (convergence check, JSON print, `forceSerial`, NLDD local solver). Justify any GPU/CPU default divergence.

## Cross-cutting passes

**Noise budget.** Cross-cutting findings apply only to lines the diff actually touches — never to surrounding code the author did not modify.

A nit has to be a real observation about a touched line. 

### C++ hygiene
- Includes: include what you use; transitive includes don't count (`<cstddef>` for `std::size_t`, `<algorithm>`, `<utility>`, `<type_traits>`, `<stdexcept>`). Remove now-unused ones. Grouping/ordering per `CONTRIBUTING.md`; a `.cpp` includes its own header first.
- `std::size_t` not bare `size_t`; `static_cast` not C-casts; safe conversion at `int`-taking boundaries (MPI, cuBLAS).
- Always initialise `bool` and scalar members. `explicit` on single-argument constructors. `override` without redundant `virtual`.
- `enum class` over `int`/`bool` flags; explicit enumerator values; no naked `bool` parameters in client-facing APIs. No `default:` in a switch over an enum — enumerate all and throw after, so the compiler still warns on new enumerators.
- `.empty()` over `.size() == 0`; `contains()` over `find() != end()`; `std::clamp`, `std::any_of`, `std::accumulate`, `std::ranges::` algorithms over hand-rolled loops. `numeric_limits<double>::lowest()`, not `min()`, for "most negative".
- `const auto&` in range-for; heavy arguments by `const&`; `std::string_view` or `const std::string&`, never `const std::string` by value.
- No raw `new`/`delete`; `std::make_unique`. `unique_ptr` over `shared_ptr` for unique ownership.
- `std::format` over concatenation and `std::to_string`. Never `std::cerr` / `printf` / `exit()` — use `OpmLog` and exceptions.
- `assert` is for programmer error only, never for user or deck input, and never as the sole guard.
- Don't reimplement Dune or OPM utilities: `Dune::DynamicMatrix`, `RegulaFalsi` in `RootFinders.hpp`, `SparseTable<int>` for ragged per-cell data, `PropertyTree` instead of raw Boost property_tree.
- Naming must not lie: `has*()` must be `const` and non-mutating; `add()` that overwrites is `assign()`; the file name matches the entity in it. Members use trailing underscore in opm-simulators, `m_` in gpuistl.

### Headers, build, CMake
- Never include a private/implementation header from an installed public header — forward-declare.
- Template pattern: declarations in `.hpp`, definitions in `_impl.hpp`, and a `.cpp` that is the *only* includer of `_impl.hpp` plus explicit instantiations. A `.cpp` without that split buys nothing.
- Every new source/header/test goes in `CMakeLists_files.cmake` in the correct list (`MAIN_SOURCE_FILES`, `PUBLIC_HEADER_FILES`, `TEST_SOURCE_FILES`, `TEST_DATA_FILES`), alphabetically; `.hpp` never in `MAIN_SOURCE_FILES`; delete entries for removed files. GPU files use `ADD_CUDA_OR_HIP_FILE`. This one is cheap to verify — grep the file list for every added filename.
- New public dependencies route through the `*_prereqs_hook`. Link imported targets (`MPI::MPI_C`, `METIS::METIS`, `fmt::fmt`, `Dune::Common`), never raw `${FOO_LIBRARIES}`.
- New files need the GPL-3.0-or-later header; moved code keeps its original copyright.

### Performance (hot paths)
- No allocation, file I/O, model loading, string construction or `std::map`/`std::set` lookup in per-cell / per-Newton-iteration loops.
- No logging inside per-cell or per-iteration paths — at most once per timestep, and only from rank 0. Unconditional `OpmLog::debug("..." + std::to_string(x))` builds the string regardless of log level; put traces behind a compile-time macro.
- `reserve()` / size up-front instead of repeated `push_back`; `emplace` in place; avoid intermediate containers; watch copies from getters (`const auto&`, not `const auto`).
- Watch per-connection / per-cell memory growth: an added `optional<>` or array is paid by every object in the model.
- `OPM_TIMEBLOCK` on genuinely hot regions only; use the `_LOCAL` variants behind `DETAILED_PROFILING` for inner loops.
- Claims of speedup or "no regression" need a number and a deck name from the author or CI — do not run benchmarks yourself.

### Tests
- New non-trivial logic, comparison/ordering predicates, and new parser behaviour all need unit tests (`tests/parser/*Tests.cpp`, `BOOST_CHECK*`).
- A test that writes a file must read it back and check the contents. A test must fail on the old code.
- Behaviour changes need a Jenkins run with `failure_report`, and **each** regression difference explained as expected, an improvement, or needing investigation — never waved through.
- Reference-solution updates need a paired opm-tests PR, referenced so both merge together.
- Prefer extending an existing fast deck over adding a slow case; stabilise a flaky case (TUNING, `--solver-max-time-step-in-days`) rather than accepting noise.

### Style and diff hygiene
- No wholesale reformatting, no whitespace-only edits to untouched lines, no files touched without a real change; it makes `git blame` and the file's history less useful.
- Rebase onto master rather than merging master into the branch. Squash conflict-marker and "minor"/"cleanup" commits; keep commits on topic with messages that explain the change.

## CI and merge protocol

Jenkins is triggered by a PR comment, not automatically:

```
jenkins build this please
jenkins build this failure_report please          # behaviour changes: get the regression PDF
jenkins build this update_data please             # reference solutions need regenerating
jenkins build this serial rocm hipify please      # GPU changes: AMD + HIP conversion
jenkins build this opm-simulators=7249 opm-tests=1020 please   # cross-repo pinning
```

Other flags seen in use: `wheels`, `ignore_extra`, `only_summary`, `no downstreams`, `clang`, `nompi`, `debug`, `hypre`, `shared`, `gcc_lto`. Pin keys: `opm-common`, `opm-simulators`, `opm-grid`, `opm-models`, `opm-upscaling`, `opm-tests`. Benchmarks are requested with `benchmark please`.

A PR merges when it is approved **and** the build check is green; where reference data changed, it must be installed on the CI system first. Cross-repo work waits for the upstream PR to merge first. Merging is a maintainer action.

## Do not

- Invent findings. If the diff doesn't show it, ask instead of asserting — and if you cannot describe a concrete path to failure, it is a Question, not a defect.
- Apply the whole checklist to a one-line fix, or comment on code the diff does not touch.
- Demand refactors unrelated to the PR's scope, or block on preferences not tied to correctness, compatibility, determinism or performance.
- Request changes over nits, or produce a wall of low-value comments.
- Present a contested style preference as a project rule. Braced initialisers versus explicit constructors, and brace placement, are not settled across the project — leave them alone.
