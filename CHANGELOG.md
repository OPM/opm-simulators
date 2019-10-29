# Changelog
All notable changes to this module (and sometimes relevant changes to other modules) will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
but this project does not use [Semantic Versioning](http://semver.org/spec/v2.0.0.html).
Instead, version numbers indicate the year and month of release.

## [2019-10] - 2019-10-31

### Added
- Support for additional keywords: ROCKCOMP keyword (water induced rock
  compaction). COMPDAT, UDQ, UDA, GRUP, FILLEPS
- Support for GOR checking for WECON
- Improved support for WTEST
- Support for one phase simulations
- Experimental support for driftCompensation (off by default, enable
  using --ecl-enable-drift-compensation=true)
- Experimental support for foam module (keywords FOAM; FOAMADS,
  FOAMFSC, FOAMMOB, FOAMOPTS, FOAMROCK, WFOAM)
- Many unsupported keywords are now listed as missing features (those
  starting with A - M, R, T, V, W, Z).
- Write well potentials to restart files if needed.
- Support different edge weights when loadbalancing (uniform,
  transmissibilities, log of transmissibilities)
- Improved support for multisegment wells (now enabled by default).
- Many new regression tests.
- Added an experimental  linear solver subsystem (including complex
  CPR solvers) that is configurable during runtime (not parallel,
  needs DUNE>=2.6)
- Use refactored well implementation from opm-common.
- Output NNCs in an Eclipse compliant manner.
- Output several diagnostics when parsing a deck.

### Fixed
- Made dune-fem version information available
- restart values are only read once (twice before)
- Fixed several bugs concerning restart
- ebos now logs to *.PRT and *DBG files
- Abort runs without reading the deck if command line parameters are incorrect.
- Use grid region mapping from opm-grid.
- Fixed negative thp values form extrapolation using VFP tables.
- Logging from well testing is now in *.PRT and *.LOG files
- Serveral fixes to multisegement well model.
- *.INIT and *.GRID files are also output on restart.
- Do not update RESV for prediction produces.
- Output FPRP insteas of ovewritung FPR values.
- Always write transmissibilities between vertical neighbours to TRANZ
  (even for NNCs).
- Support upcoming DUNE 2.7

### Renaming / Restructuring
- former module ewoms was renamed to opm-models. With this move
  several directories have been renamed (ewoms/blackoil ->
  opm/models/blackoil, ewoms/common ->  opm/models/utils,
  ewoms/disc -> opm/models/discretization, ewoms/io -> opm/models/io,
  ewoms->parallel -> opm/models/parallel, ewoms/linear ->
  opm/simulators/linalg). The namespace Ewoms was renamed to Opm.
- Most files from opm/autodiff/ have been moved to
  opm/simulators/{aquifers,linalg,utils,wells}
  
### Removed
- Dependency on libecl

## [2018.10 - 2019.04] is skipped.

## [2018.04] - 2018-04.30

### Added
- More integration tests for Flow, in particular for multisegment wells, solvent and polymer.
- Support placing the TUNING keyword (first line only, time limits for simulator steps) in the SCHEDULE section.
- Support RESV-controlled wells for solvent runs.
- More output options supported, such as dissolved gas and reservoir volume rates.
- Support two-phase oil-water with polymer without a dummy gas phase.
- Experimental CPR preconditioner is available for Flow.

### Changed
- Some parts of the former opm-core module has been moved to this module.
- Refactoring: well models are now more independent and self-contained.
- Refactoring: the Schedule objects are now independent of the EclipseState objects and can be modified from the initial deck-given state.
- Refactoring: EQUIL initialization has been refactored and moved to the ewoms module.
- Refactoring: use refactored output facilities.
- Refactoring: more flexible treatment of linear solver.
- Improved and cleaned up log and PRT output from Flow.
- Changed some default tuning parameters for Flow: default initial timestep size is 1 day, maximum pressure change per iteration is limited, allow setting the relaxed convergence tolerance.
- Other minor refactoring: SimulatorBlackoilEbos::run(), BlackoilModelEbos::prepareStep() and afterStep(), SimulatorReport usage, cleanup component and phase index usage.
- Allow using Dune 2.6.

### Fixed
- Fix bugs related to running Flow in parallel that caused slightly wrong results and bad performance in some cases, in particular the "model 2" case with 8 threads.
- Fix bugs for polymer simulation, in particular relating to partially mixed polymer and shear-dependent behaviour.
- Fix a bug related to wells with a zero rate target.
- Fixed various minor bugs, including bugs in equilibration initialization, printing very long timesteps, single-phase RESV injection controls, RESV controls in the presence of undersaturated fluids, logging for old two-phase simulators.

### Removed
- The stand-alone ebos-based (i.e. using the new assembly approach) Flow variants as well as most of the legacy variants have been removed, use Flow instead.
- The legacy Flow variants for polymer, solvent etc. have been removed, since that functionality is all present in the main Flow executable.
- The old implicit solver for incompressible two-phase transport has been removed, the ad-based version is simply better.

## [2017.10] - 2017-10-31

### Added
- Added much-improved support for multisegment wells to Flow.
- Two-phase gas-oil cases are now supported in addition to water-oil (gas-water is still not supported).
- Added an experimental sequential splitting simulator using a reordering transport solver.

### Changed
- Major refactoring of well models.
- All Flow variant models have better performance. In particular this is significant for the two-phase, solvent and polymer variants, but also the basic black-oil model is faster.
- All Flow variants have now been collected in a single executable. The old, separate executables are still there, but will be removed.
- Improved log messages and reporting.
- Changed convergence testing for improved speed and stability control.
- Improved integration test coverage.
- Minor refactoring of several classes and systems.
- Build system fixes and cleanup.
- Improved Jenkins integration for data updating and static code analysis.
- Require Dune version at least 2.4, allow use of Dune 2.5.

### Fixed
- Fixed various minor bugs, including bugs in RESV controls, group controls, use of well efficiency factors, VFP treatment, fluid-in-place output.



## [2017.04] - 2017-04-26

### Added
- A Docker container has been created on Docker Hub, for easy deployment.
- Many more tests are now being run automatically to ensure correct simulator behaviour.

### Changed
- Flow now uses a new code for assembly of residuals and Jacobian matrices, resulting in significant performance gains. The new code is still based on automatic differentiation, but uses a more local approach than the previous method.
- Well handling has seen significant improvement, especially group controls. Flow is now able to correctly change well control depending on group limits in many cases previously not handled, although Flow still does not handle group control sets spanning more than two levels well. THP limits and VFP tables are handled in a more robust manner, and LRAT controls now work.
- Log output to PRT file and terminal is more precise, comprehensive and user friendly, with improved log messages and classification of messages. Timing and iteration reports have been cleaned up and now also include iterations that were wasted due to convergence failures. Also some minor output improvements for the 2-phase case have been made.
- Output possibilities for summary and restart files are now even more comprehensive, for example bubble point pressures can now be written. Restarting Flow simulations from restart files is improved, and reproduces initial runs much better than before.
- All grid-related parts of opm-core have been moved to the opm-grid module.
- Grids with inverted Z coordinates and left-handed coordinate system are now handled correctly.

### Fixed
- Many small bugs have been fixed. In total more than two hundred pull requests have been merged since last release.








