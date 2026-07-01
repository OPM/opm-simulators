# The TypeTag / Property System in opm-simulators

This document gives a broad overview of how the **TypeTag / property
system** is used in `opm-simulators` (and the underlying `opm-models`
framework). It covers:

1. What a TypeTag is and how the system works.
2. How to add a new TypeTag.
3. How to add a new property tag.
4. How to set a value (specialise a property) for a TypeTag.
5. Where the central definitions live and what to include.
6. A specification of the **blackoil (default `FlowProblem`)** TypeTag.
7. A specification of the **CO2STORE** configuration (which uses the
   blackoil TypeTag with runtime customisation rather than a separate
   TypeTag).

> The implementation lives in
> `opm-simulators/opm/models/utils/propertysystem.hh`. The handbook
> chapter (`opm-simulators/doc/handbook/propertysystem.tex`) still uses
> the legacy `NEW_TYPE_TAG`/`NEW_PROP_TAG` macro syntax; the current
> code uses C++17 struct-based syntax. **The struct-based form is the
> only one used in the codebase today.**

---

## 1. Concepts

The TypeTag system is a compile-time *traits-with-inheritance* system:

* A **TypeTag** is an empty struct in `Opm::Properties::TTag` that
  optionally lists a `using InheritsFrom = std::tuple<...>;` of parent
  TypeTags. It is the "node" in the inheritance hierarchy.
* A **Property** is a class template `template<class TypeTag, class
  MyTypeTag> struct Foo { ... };` in `Opm::Properties`. The first
  parameter is the *concrete* TypeTag the user is asking about; the
  second is the TypeTag at which the property is currently being
  searched for during the recursive lookup.
* The **default declaration** of a property uses
  `using type = UndefinedProperty;` so that the lookup machinery can
  detect whether a value has been set.
* A **value** for a property on a TypeTag is provided by a
  **partial specialisation** of the property template for that TypeTag.
* Values are queried with the helpers in `propertysystem.hh`:
  * `GetPropType<TypeTag, Properties::Foo>` — the `::type` member
    (replaces the old `GET_PROP_TYPE` macro).
  * `GetProp<TypeTag, Properties::Foo>` — the whole property struct
    (replaces the old `GET_PROP` macro).
  * `getPropValue<TypeTag, Properties::Foo>()` — the `::value` member
    for value-properties (replaces `GET_PROP_VALUE`).
* **Splices** (`Opm::Properties::Splices<TypeTag, MyTypeTag>`) allow a
  TypeTag to dynamically choose a parent at compile time (used e.g. to
  switch between `EcfvDiscretization` and `VcfvDiscretization`).

The lookup walks the `InheritsFrom` graph in depth-first order until
it finds a defined specialisation. Because the recursion is on the
*hierarchy*, a derived TypeTag can override any property from any
ancestor.

There are **no macros** anymore (`NEW_TYPE_TAG`, `NEW_PROP_TAG`,
`SET_TYPE_PROP`, `SET_INT_PROP`, …): everything is plain C++17 class
templates and `std::tuple`.

---

## 2. How to add a new TypeTag

A TypeTag is a struct in the namespace `Opm::Properties::TTag`. It
must be visible everywhere that properties are specialised for it.

* **Where to put it:** Either inline in the `.cpp`/`.hpp` that needs
  it (see `flow/flow_gaswater_dissolution.cpp`, lines ~37–40), or in a
  small dedicated header so that several files can share it (see
  `opm/simulators/flow/TTagFlowProblemTPFA.hpp`,
  `TTagFlowProblemGasWater.hpp`, `TTagFlowProblemOnePhase.hpp`).
* **What to include:** `<tuple>` plus a forward declaration of any
  parent TypeTag struct used in `InheritsFrom`. The full definition of
  a parent only needs to be visible at the point where you *specialise
  a property* against it.
* **Pattern:**

  ```cpp
  // namespace Opm::Properties::TTag
  struct ParentTag;                  // forward declaration is enough
  struct MyNewTag {
      using InheritsFrom = std::tuple<ParentTag /*, OtherParent ...*/>;
  };
  ```

  Multiple inheritance is allowed and ordered: earlier entries in the
  tuple win when the same property is set in two ancestors.
* **Naming:** Existing TypeTags are CamelCase, with `FlowProblem...`
  for the simulator-level ones, and `BlackOilModel`, `MultiPhaseBaseModel`,
  `EcfvDiscretization`, `VcfvDiscretization`, `NumericModel`, … for
  framework-level ones.
* **Forward declaration trick:** Inside header files in
  `opm/simulators/flow/TTag*.hpp`, the parent `FlowProblem` is only
  forward declared. The full `FlowProblem` TypeTag is defined in
  `opm/simulators/flow/BlackoilModelProperties.hpp`. This keeps the
  `TTag*` headers cheap to include.

---

## 3. How to add a new property tag

A property tag is a class template with two template parameters and a
default `UndefinedProperty` body. It must live in `Opm::Properties`.

* **Where to put it:** In a `*Properties.hpp` header co-located with
  the code that consumes it, e.g.
  `opm/simulators/flow/FlowBaseProblemProperties.hpp`,
  `opm/simulators/flow/BlackoilModelProperties.hpp`,
  `opm/models/utils/basicproperties.hh`,
  `opm/models/blackoil/blackoilproperties.hh`,
  `opm/models/discretization/common/fvbaseproperties.hh`.
* **What to include:** `<opm/models/utils/propertysystem.hh>`.
* **Pattern (type-valued property):**

  ```cpp
  namespace Opm::Properties {
  template<class TypeTag, class MyTypeTag>
  struct AquiferModel { using type = UndefinedProperty; };
  }
  ```

* **Pattern (value-valued property):** the convention is the same; the
  *specialisations* later expose either `using type = ...;` or
  `static constexpr T value = ...;`. The default body should still use
  `UndefinedProperty` so that the lookup detects "unset". For
  enum-valued properties (e.g. `EnergyModuleType`), use
  `using type = Properties::UndefinedProperty;` in the default and
  `static constexpr Foo value = ...;` in the specialisations.
* **Recommendation:** A property tag may be re-declared in any file
  that uses it, but the de-facto convention in this code base is to
  declare each tag once in the most relevant `*Properties.hpp` header
  and `#include` that header.

---

## 4. How to set a value for a property on a TypeTag

You provide a **partial specialisation** of the property template for
`(TypeTag, MyTypeTag) = (TypeTag, TTag::YourTag)`. The first parameter
is left as a template parameter so the property still depends on the
concrete TypeTag for transitive lookups.

* **Type-valued:**

  ```cpp
  template<class TypeTag>
  struct Problem<TypeTag, TTag::FlowProblem>
  { using type = FlowProblemBlackoil<TypeTag>; };
  ```

* **bool-valued:**

  ```cpp
  template<class TypeTag>
  struct EnableBrine<TypeTag, TTag::FlowProblem>
  { static constexpr bool value = false; };
  ```

* **Numeric / enum-valued:**

  ```cpp
  template<class TypeTag>
  struct EnergyModuleType<TypeTag, TTag::FlowGasWaterDissolutionProblem>
  { static constexpr EnergyModules value = EnergyModules::ConstantTemperature; };
  ```

* **Computed type:** the body may use `GetPropType<TypeTag, ...>` to
  derive the value from other properties (see e.g. `Indices` for
  `BlackOilModel` in `opm/models/blackoil/blackoilmodel.hh`).
* **Where to put it:**
  * Set it next to where the *TypeTag* is defined if it is fundamental
    to that TypeTag.
  * Set it next to where the *property tag* is defined if it is a
    cross-cutting default for many TypeTags.
  * Set it inside the `.cpp` that produces an executable when the
    setting is specific to that executable (e.g. the `Linearizer` and
    `LocalResidual` overrides in `flow/flow_blackoil.cpp`).

### Querying a property in code

```cpp
using Scalar     = GetPropType<TypeTag, Properties::Scalar>;
using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
constexpr bool brine = getPropValue<TypeTag, Properties::EnableBrine>();
```

---

## 5. Where the definitions live

The TypeTag/property system spans `opm-models` (framework) and
`opm-simulators` (`flow` simulator and friends). Useful entry points:

| Area | Path |
|------|------|
| Property system core | `opm-simulators/opm/models/utils/propertysystem.hh` |
| Basic property tags + `NumericModel` / `ImplicitModel` TypeTags | `opm-simulators/opm/models/utils/basicproperties.hh` |
| Discretization TypeTags (`EcfvDiscretization`, `VcfvDiscretization`) | `opm-simulators/opm/models/discretization/{ecfv,vcfv}/*properties.hh`, `opm/models/discretization/common/fvbaseproperties.hh` |
| Multiphase base TypeTag (`MultiPhaseBaseModel`) and splices | `opm-simulators/opm/models/common/multiphasebasemodel.hh` |
| Black-oil model TypeTag (`BlackOilModel`) and properties | `opm-simulators/opm/models/blackoil/blackoilmodel.hh`, `blackoilproperties.hh` |
| Other physics models | `opm-simulators/opm/models/{immiscible,ncp,pvs,flash,ptflash,richards,discretefracture}/*model.hh` |
| `flow` simulator base TypeTag (`FlowBaseProblem`) | `opm-simulators/opm/simulators/flow/FlowBaseProblemProperties.hpp` |
| `flow` blackoil base (`FlowBaseProblemBlackoil`) | `opm-simulators/opm/simulators/flow/FlowProblemBlackoilProperties.hpp` |
| Default `flow` blackoil TypeTag (`FlowProblem`) | `opm-simulators/opm/simulators/flow/BlackoilModelProperties.hpp` |
| Light-weight TypeTag headers (forward-declare `FlowProblem`, define a child) | `opm-simulators/opm/simulators/flow/TTagFlowProblem*.hpp` |
| Per-executable TypeTags & overrides | `opm-simulators/flow/flow_*.cpp` |
| Composition / experimental TypeTags | `opm-simulators/flowexperimental/*` and `opm/simulators/flow/FlowProblemCompProperties.hpp` |

### Includes you typically need

```cpp
#include <opm/models/utils/propertysystem.hh>          // GetPropType, GetProp, getPropValue
#include <opm/models/utils/basicproperties.hh>         // NumericModel, Scalar, Vanguard, ...
#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>     // FlowBaseProblem
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp> // FlowBaseProblemBlackoil
#include <opm/simulators/flow/BlackoilModelProperties.hpp>       // FlowProblem (default blackoil)
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>           // FlowProblemTPFA
#include <tuple>                                                 // std::tuple for InheritsFrom
```

For `flow`-style executables, also include the simulator main:
`<opm/simulators/flow/Main.hpp>` and
`<opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>`.

---

## 6. The blackoil (default) TypeTag specification

The default `flow` blackoil simulator (`flow/flow_blackoil.cpp`) uses
`Opm::Properties::TTag::FlowProblemTPFA` which inherits from
`FlowProblem`. The relevant chain, from most-derived to most-base, is:

```
FlowProblemTPFA            (opm/simulators/flow/TTagFlowProblemTPFA.hpp)
   └── FlowProblem         (opm/simulators/flow/BlackoilModelProperties.hpp)
        ├── FlowBaseProblemBlackoil  (opm/simulators/flow/FlowProblemBlackoilProperties.hpp)
        │     └── FlowBaseProblem    (opm/simulators/flow/FlowBaseProblemProperties.hpp)
        │          └── CpGridVanguard     (opm/simulators/flow/CpGridVanguard.hpp)
        └── BlackOilModel       (opm/models/blackoil/blackoilmodel.hh)
              ├── VtkBlackOilPolymer
              └── MultiPhaseBaseModel    (opm/models/common/multiphasebasemodel.hh)
                    └── splices: SpatialDiscretizationSplice (default EcfvDiscretization)
                           └── ImplicitModel
                                  └── NumericModel  (opm/models/utils/basicproperties.hh)
```

### TypeTag definitions (verbatim)

* `FlowProblemTPFA` — `opm/simulators/flow/TTagFlowProblemTPFA.hpp`:
  ```cpp
  struct FlowProblem;                       // forward decl
  struct FlowProblemTPFA {
      using InheritsFrom = std::tuple<FlowProblem>;
  };
  ```
* `FlowProblem` — `opm/simulators/flow/BlackoilModelProperties.hpp`:
  ```cpp
  struct FlowProblem {
      using InheritsFrom = std::tuple<FlowBaseProblemBlackoil, BlackOilModel>;
  };
  ```
* `FlowBaseProblemBlackoil` — `opm/simulators/flow/FlowProblemBlackoilProperties.hpp`:
  ```cpp
  struct FlowBaseProblemBlackoil {
      using InheritsFrom = std::tuple<FlowBaseProblem>;
  };
  ```
* `FlowBaseProblem` — `opm/simulators/flow/FlowBaseProblemProperties.hpp`:
  ```cpp
  struct FlowBaseProblem {
      using InheritsFrom = std::tuple<CpGridVanguard>;
  };
  ```

### Key property values bound to the blackoil TypeTag

Set on `TTag::FlowBaseProblemBlackoil` (in `FlowProblemBlackoilProperties.hpp`):

| Property | Value | Notes |
|---|---|---|
| `NonlinearSystem` | `BlackoilModel<TypeTag>` | |
| `Problem` | `FlowProblemBlackoil<TypeTag>` | The application problem class |
| `Model` | `FIBlackOilModel<TypeTag>` | Fully-implicit blackoil model |
| `FluxModule` | `NewTranFluxModule<TypeTag>` | Eclipse "NEWTRAN" transmissibilities |
| `MaterialLaw` (set on `FlowBaseProblem`) | `EclMaterialLaw::Manager<Traits>::MaterialLaw` | Three-phase, with hysteresis/endpoint scaling toggled by other properties |

Set on `TTag::FlowProblem` (in `BlackoilModelProperties.hpp`):

| Property | Value |
|---|---|
| `BlackoilConserveSurfaceVolume` | `true` (formulate equations in surface volumes) |
| `UseVolumetricResidual` | `false` |
| `AquiferModel` | `BlackoilAquiferModel<TypeTag>` |
| `WellModel` | `BlackoilWellModel<TypeTag>` |
| `LinearSolverSplice` | `TTag::FlowIstlSolver` |
| `EnableConvectiveMixing` | `true` |
| `EnablePolymer` | `false` |
| `EnableSolvent` | `false` |
| `EnableFoam` | `false` |
| `EnableBrine` | `false` |
| `EnableSaltPrecipitation` | `false` |
| `EnableBioeffects` | `false` |
| `EnableDispersion` | `false` |
| `EnableDebuggingChecks` | `false` |

Set on `TTag::FlowProblemTPFA` (overrides applied only for the
TPFA-assembled default executable, in `flow/flow_blackoil.cpp`):

| Property | Value |
|---|---|
| `Linearizer` | `TpfaLinearizer<TypeTag>` |
| `LocalResidual` | `BlackOilLocalResidualTPFA<TypeTag>` |
| `EnableDiffusion` | `false` |
| `AvoidElementContext` | `true` |

Inherited from `TTag::BlackOilModel`
(`opm/models/blackoil/blackoilmodel.hh`), the most important values are:

| Property | Value |
|---|---|
| `LocalResidual` | `BlackOilLocalResidual<TypeTag>` (overridden by `FlowProblemTPFA`) |
| `NewtonMethod` | `BlackOilNewtonMethod<TypeTag>` |
| `Model` | `BlackOilModel<TypeTag>` (overridden by `FlowBaseProblemBlackoil` to `FIBlackOilModel<TypeTag>`) |
| `BaseProblem` | `BlackOilProblem<TypeTag>` |
| `RateVector` | `BlackOilRateVector<TypeTag>` |
| `BoundaryRateVector` | `BlackOilBoundaryRateVector<TypeTag>` |
| `PrimaryVariables` | `BlackOilPrimaryVariables<TypeTag>` |
| `IntensiveQuantities` | `BlackOilIntensiveQuantities<TypeTag>` |
| `ExtensiveQuantities` | `BlackOilExtensiveQuantities<TypeTag>` |
| `FluxModule` | `BlackOilDarcyFluxModule<TypeTag>` (overridden by blackoil base to `NewTranFluxModule<TypeTag>`) |
| `Indices` | `BlackOilVariableAndEquationIndices<EnableSolvent, EnableExtbo, EnablePolymer, EnergyAsBool, EnableFoam, EnableBrine, /*PVOffset*/ 0, EnableBioeffects>` |
| `FluidSystem` | `BlackOilFluidSystem<Scalar>` |
| Module enables (`EnableSolvent`, `EnableExtbo`, `EnablePolymer`, `EnablePolymerMW`, `EnableFoam`, `EnableBrine`, `EnableVapwat`, `EnableDisgasInWater`, `EnableSaltPrecipitation`, …) | all `false` by default |

Inherited further from `TTag::MultiPhaseBaseModel`: `NumPhases`,
`NumEq`, `NumDerivatives`, the `SpatialDiscretizationSplice` (defaults
to `EcfvDiscretization`), gravity-related defaults, etc.

The grid comes from `TTag::CpGridVanguard` via `FlowBaseProblem`,
giving Eclipse-style corner-point grids by default.

### Variants (executables) built on top of `FlowProblem`

Each `flow/flow_*.cpp` declares its own TypeTag inheriting from
`FlowProblem` and flips a small set of properties. Typical examples:

* `flow_blackoil.cpp` → `FlowProblemTPFA` (TPFA assembly + default
  blackoil).
* `flow_gaswater_dissolution.cpp` → `FlowGasWaterDissolutionProblem`,
  enabling `EnableDisgasInWater` and `EnableVapwat`, switching `Indices`
  to `BlackOilTwoPhaseIndices`, etc.
* `flow_brine.cpp`, `flow_polymer.cpp`, `flow_solvent.cpp`,
  `flow_foam.cpp`, `flow_brine_saltprecipitation.cpp`,
  `flow_blackoil_temp.cpp`, etc. — each sets the corresponding
  `Enable*` property to `true` and may swap `Indices`/`FluidSystem`.
* `flow_blackoil_alugrid.cpp`, `flow_blackoil_polyhedralgrid.cpp` —
  swap the `Vanguard`/`Grid` properties.

---

## 7. The CO2STORE "TypeTag" specification

> **Important:** there is **no dedicated `CO2Store` TypeTag** in
> `opm-simulators`. CO2 storage runs use the **same `FlowProblem` /
> `FlowProblemTPFA` TypeTag as the default `flow` blackoil simulator**;
> the CO2-storage behaviour is selected *at runtime* via the Eclipse
> deck keyword `CO2STORE` and queried from the runspec
> (`eclState.runspec().co2Storage()`). The same is true for `H2STORE`.

This is a deliberate design choice: the blackoil model is reinterpreted
so that the "oil" component plays the role of brine/water and the "gas"
component plays the role of CO₂; thereby the existing fully-implicit
black-oil machinery is reused.

### How it is selected

* **At deck level:** add `CO2STORE` (and optionally `THERMAL`,
  `DISGASW`, `VAPWAT`, …) in the `RUNSPEC` section. The standard
  `flow` binary will then route the simulation through the CO₂ storage
  code paths.
* **At simulator level:** the same TypeTag (`FlowProblemTPFA` for the
  default `flow` executable, or `FlowGasWaterDissolutionProblem`,
  `FlowGasWaterEnergyProblem`, … for the dedicated executables) is
  used; the chosen TypeTag mainly controls *which Indices and which
  Enable* flags are compile-time true. CO₂STORE-specific behaviour is
  then activated dynamically.

### TypeTags typically used to run CO2STORE decks

| Deck features | Recommended `flow` executable | TypeTag |
|---|---|---|
| Plain 3-phase blackoil with `CO2STORE` | `flow` | `FlowProblemTPFA` (inherits `FlowProblem`) |
| Gas–water with dissolution / vaporisation | `flow_gaswater_dissolution` | `FlowGasWaterDissolutionProblem` (inherits `FlowProblem`, sets `EnableDisgasInWater = true`, `EnableVapwat = true`, two-phase `Indices`) |
| Gas–water with energy | `flow_gaswater_energy` | `FlowGasWaterEnergyProblem` (inherits `FlowProblem`, enables thermal energy module) |
| Gas–water with brine / salt precipitation | `flow_gaswater_brine`, `flow_gaswater_saltprec_*` | corresponding `FlowGasWater*Problem` TypeTag, enabling `EnableBrine` / `EnableSaltPrecipitation` |
| Brine + salt precipitation (3-phase) | `flow_brine_saltprecipitation` | `FlowBrineSaltPrecipitationProblem` |
| With diffusion | any of the above with `EnableDiffusion = true` set on the TypeTag (see `flow_gaswater_dissolution_diffuse.cpp`) | dedicated `*Diffuse` TypeTag |

In every case the TypeTag is one of the children of `FlowProblem` and
therefore inherits the entire blackoil property stack documented in
section 6.

### Properties that switch behaviour for CO2 storage

The following compile-time properties are usually flipped on the
TypeTag chosen for a CO2STORE simulation:

| Property | Typical value for CO2STORE | Set in |
|---|---|---|
| `EnableDisgasInWater` | `true` (CO₂ dissolves in water) | child TypeTag, e.g. `FlowGasWaterDissolutionProblem` |
| `EnableVapwat` | `true` (water vaporises into the gas/CO₂ phase) | child TypeTag |
| `EnableDiffusion` | `true` for `*_diffuse` variants, `false` otherwise | child TypeTag |
| `EnableBrine` | `true` if salt is modelled | child TypeTag |
| `EnableSaltPrecipitation` | `true` for salt-precipitation runs | child TypeTag |
| `EnergyModuleType` | `EnergyModules::FullyImplicitThermal` for energy/thermal CO2STORE; `EnergyModules::ConstantTemperature` otherwise | child TypeTag |
| `Indices` | `BlackOilTwoPhaseIndices<...>` for gas-water runs (oil component disabled via `disabledCompIdx = FluidSystem::oilCompIdx`) | child TypeTag |

The runtime switches consumed inside the model are the deck-controlled
flags, e.g.:

* `eclState.runspec().co2Storage()` — picks the CO2-storage variant of
  the blackoil fluid system, sets a default isothermal temperature for
  isothermal runs (see `FlowProblemBlackoil.hpp`, around the
  `// For CO2STORE and H2STORE` comment), reuses gas-phase storage
  formulas for CO₂ (see `blackoilintensivequantities.hh`,
  `OutputBlackoilModule.hpp`, `BlackoilWellModel.hpp`,
  `AquiferInterface.hpp`, `RelpermDiagnostics.cpp`,
  `InitStateEquil_impl.hpp`, `MultisegmentWell_impl.hpp`,
  `StandardWell_impl.hpp`, …).
* `eclState.runspec().h2Storage()` — analogous for H₂.

### Concrete property "specification" of a typical CO2STORE setup

Taking `flow` (default executable) running a `CO2STORE` deck, the
**effective TypeTag** is `Opm::Properties::TTag::FlowProblemTPFA`
(defined in
`opm/simulators/flow/TTagFlowProblemTPFA.hpp`). Its compile-time
property values are exactly those documented in section 6, i.e.:

* `Problem = FlowProblemBlackoil<TypeTag>`,
  `Model = FIBlackOilModel<TypeTag>`,
* `LocalResidual = BlackOilLocalResidualTPFA<TypeTag>`,
* `Linearizer = TpfaLinearizer<TypeTag>`,
* `FluidSystem = BlackOilFluidSystem<Scalar>` (the CO₂-storage
  *parameters* of this fluid system are configured at runtime from the
  deck — the *type* is the same),
* `WellModel = BlackoilWellModel<TypeTag>`,
* `AquiferModel = BlackoilAquiferModel<TypeTag>`,
* All `Enable*` extension flags `false`,
* `Indices = BlackOilVariableAndEquationIndices<…all false…, /*PVOffset*/0, /*EnableBioeffects*/false>`.

Running the same deck through, e.g., `flow_gaswater_dissolution`
changes the TypeTag to `FlowGasWaterDissolutionProblem`, which on top
sets `EnableDisgasInWater = true`, `EnableVapwat = true`,
`EnergyModuleType = ConstantTemperature`, and replaces `Indices` with a
two-phase variant where the oil component is disabled.

### Summary

* "blackoil (default)" maps to **one and the same TypeTag** as the
  `flow` binary uses: `TTag::FlowProblemTPFA`, inheriting from
  `TTag::FlowProblem` and ultimately the full blackoil/multiphase
  property stack.
* "co2store" is **not a separate TypeTag** but a *configuration* of
  that same TypeTag (or one of its gas-water siblings) selected by the
  deck keyword `CO2STORE`. The TypeTag mostly controls which
  compile-time `Enable*` flags and which `Indices` are used; the
  CO₂-vs-blackoil distinction is made at runtime through
  `eclState.runspec().co2Storage()`.

---

## Quick-reference checklist for new TypeTags / properties

1. **New TypeTag for a new executable / variant:**
   * Pick the closest existing TypeTag (almost always `FlowProblem` for
     `flow`-style executables).
   * Add `struct MyTag { using InheritsFrom = std::tuple<FlowProblem>; };`
     in the `.cpp` (or in a `TTag*.hpp` header if shared).
   * Add the `Linearizer`, `LocalResidual`, `Indices`, `Enable*`, …
     specialisations needed for your physics.
   * In the `main`-style entry point use
     `FlowMain<Properties::TTag::MyTag>` /
     `mainObject->runStatic<Properties::TTag::MyTag>()`.

2. **New property tag:**
   * Declare `template<class TypeTag, class MyTypeTag> struct Foo
     { using type = UndefinedProperty; };` in the most relevant
     `*Properties.hpp` header under `namespace Opm::Properties`.
   * Provide a default value on the most-base TypeTag where it makes
     sense (e.g. on `BlackOilModel` or `FlowProblem`).
   * Consume it via `GetPropType<TypeTag, Properties::Foo>` /
     `getPropValue<TypeTag, Properties::Foo>()`.

3. **New value for an existing property on an existing TypeTag:**
   * Partial-specialise the property template:
     `template<class TypeTag> struct Foo<TypeTag, TTag::Bar> { ... };`
     in a header included by every translation unit that resolves
     `Foo` against `TTag::Bar` (typically next to the `TTag::Bar`
     definition itself).

