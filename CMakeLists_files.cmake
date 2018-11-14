# This file sets up five lists:
# MAIN_SOURCE_FILES     List of compilation units which will be included in
#                       the library. If it isn't on this list, it won't be
#                       part of the library. Please try to keep it sorted to
#                       maintain sanity.
#
# TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
# TEST_DATA_FILES       Files from the source three that should be made
#                       available in the corresponding location in the build
#                       tree in order to run tests there.
#
# EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#                       build, but which is not part of the library nor is
#                       run as tests.
#
# PUBLIC_HEADER_FILES   List of public header files that should be
#                       distributed together with the library. The source
#                       files can of course include other files than these;
#                       you should only add to this list if the *user* of
#                       the library needs it.

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
  opm/autodiff/Compat.cpp
  opm/autodiff/ExtractParallelGridInformationToISTL.cpp
  opm/autodiff/NewtonIterationBlackoilInterleaved.cpp
  opm/autodiff/NewtonIterationUtilities.cpp
  opm/autodiff/GridHelpers.cpp
  opm/autodiff/moduleVersion.cpp
  opm/autodiff/BlackoilPropsAdFromDeck.cpp
  opm/autodiff/BlackoilModelParameters.cpp
  opm/autodiff/WellDensitySegmented.cpp
  opm/autodiff/LinearisedBlackoilResidual.cpp
  opm/autodiff/MPIUtilities.cpp
  opm/autodiff/VFPProdProperties.cpp
  opm/autodiff/VFPProdPropertiesLegacy.cpp
  opm/autodiff/VFPInjProperties.cpp
  opm/autodiff/VFPInjPropertiesLegacy.cpp
  opm/autodiff/MissingFeatures.cpp
  opm/core/props/rock/RockFromDeck.cpp
  opm/core/props/satfunc/RelpermDiagnostics.cpp
  opm/core/props/satfunc/SaturationPropsFromDeck.cpp
  opm/core/simulator/BlackoilState.cpp
  opm/core/simulator/SimulatorReport.cpp
  opm/core/utility/Event.cpp
  opm/core/wells/InjectionSpecification.cpp
  opm/core/wells/ProductionSpecification.cpp
  opm/core/wells/WellCollection.cpp
  opm/core/wells/WellsGroup.cpp
  opm/core/wells/WellsManager.cpp
  opm/core/wells/well_controls.c
  opm/core/wells/wells.c
  opm/polymer/PolymerBlackoilState.cpp
  opm/simulators/WellSwitchingLogger.cpp
  opm/simulators/timestepping/TimeStepControl.cpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.cpp
  opm/simulators/timestepping/SimulatorTimer.cpp
  opm/simulators/timestepping/gatherConvergenceReport.cpp
  )

if(PETSc_FOUND)
  list(APPEND MAIN_SOURCE_FILES opm/core/linalg/LinearSolverPetsc.cpp)
endif()


# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/test_autodiffhelpers.cpp
  tests/test_autodiffmatrix.cpp
  tests/test_blackoil_amg.cpp
  tests/test_block.cpp
  tests/test_convergencereport.cpp
  tests/test_graphcoloring.cpp
  tests/test_span.cpp
  tests/test_syntax.cpp
  tests/test_scalar_mult.cpp
  tests/test_vfpproperties.cpp
  tests/test_milu.cpp
  tests/test_multmatrixtransposed.cpp
  tests/test_wellmodel.cpp
#  tests/test_thresholdpressure.cpp
  tests/test_wellswitchlogger.cpp
  tests/test_timer.cpp
  tests/test_invert.cpp
  tests/test_event.cpp
  tests/test_wells.cpp
  tests/test_equil_legacy.cpp
  tests/test_blackoilstate.cpp
  tests/test_wellsmanager.cpp
  tests/test_wellcontrols.cpp
  tests/test_wellsgroup.cpp
  tests/test_wellcollection.cpp
  tests/test_stoppedwells.cpp
  tests/test_relpermdiagnostics.cpp
  tests/test_norne_pvt.cpp
  )

if(MPI_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_parallelistlinformation.cpp)
endif()

list (APPEND TEST_DATA_FILES
  tests/VFPPROD1
  tests/VFPPROD2
  tests/msw.data
  tests/TESTTIMER.DATA
  tests/TESTWELLMODEL.DATA
  tests/liveoil.DATA
  tests/capillary.DATA
  tests/capillary_overlap.DATA
  tests/capillarySwatinit.DATA
  tests/deadfluids.DATA
  tests/equil_livegas.DATA
  tests/equil_liveoil.DATA
  tests/equil_rsvd_and_rvvd.DATA
  tests/wetgas.DATA
  tests/satfuncEPS_B.DATA
  tests/testBlackoilState1.DATA
  tests/testBlackoilState2.DATA
  tests/wells_manager_data.data
  tests/wells_manager_data_expanded.data
  tests/wells_manager_data_wellSTOP.data
  tests/wells_group.data
  tests/wells_stopped.data
  tests/relpermDiagnostics.DATA
  tests/norne_pvt.data
  )


# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
  examples/find_zero.cpp
  )

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/autodiff/AutoDiffBlock.hpp
  opm/autodiff/AutoDiffHelpers.hpp
  opm/autodiff/AutoDiffMatrix.hpp
  opm/autodiff/AutoDiff.hpp
  opm/autodiff/BlackoilAmg.hpp
  opm/autodiff/BlackoilDetails.hpp
  opm/autodiff/BlackoilLegacyDetails.hpp
  opm/autodiff/BlackoilModel.hpp
  opm/autodiff/BlackoilModelBase.hpp
  opm/autodiff/BlackoilModelBase_impl.hpp
  opm/autodiff/BlackoilModelEnums.hpp
  opm/autodiff/BlackoilModelParameters.hpp
  opm/autodiff/BlackoilModelParametersEbos.hpp
  opm/autodiff/BlackoilPressureModel.hpp
  opm/autodiff/BlackoilPropsAdFromDeck.hpp
  opm/autodiff/Compat.hpp
  opm/autodiff/CPRPreconditioner.hpp
  opm/autodiff/createGlobalCellArray.hpp
  opm/autodiff/DefaultBlackoilSolutionState.hpp
  opm/autodiff/BlackoilSequentialModel.hpp
  opm/autodiff/BlackoilReorderingTransportModel.hpp
  opm/autodiff/BlackoilTransportModel.hpp
  opm/autodiff/fastSparseOperations.hpp
  opm/autodiff/DebugTimeReport.hpp
  opm/autodiff/DuneMatrix.hpp
  opm/autodiff/ExtractParallelGridInformationToISTL.hpp
  opm/autodiff/FlowLinearSolverParameters.hpp
  opm/autodiff/FlowMain.hpp
  opm/autodiff/FlowMainEbos.hpp
  opm/autodiff/FlowMainSequential.hpp
  opm/autodiff/GeoProps.hpp
  opm/autodiff/GraphColoring.hpp
  opm/autodiff/GridHelpers.hpp
  opm/autodiff/GridInit.hpp
  opm/autodiff/ISTLSolver.hpp
  opm/autodiff/ISTLSolverEbos.hpp
  opm/autodiff/IterationReport.hpp
  opm/autodiff/moduleVersion.hpp
  opm/autodiff/NewtonIterationBlackoilInterface.hpp
  opm/autodiff/NewtonIterationBlackoilInterleaved.hpp
  opm/autodiff/NewtonIterationUtilities.hpp
  opm/autodiff/NonlinearSolver.hpp
  opm/autodiff/NonlinearSolver_impl.hpp
  opm/autodiff/NonlinearSolverEbos.hpp
  opm/autodiff/LinearisedBlackoilResidual.hpp
  opm/autodiff/ParallelDebugOutput.hpp
  opm/autodiff/ParallelOverlappingILU0.hpp
  opm/autodiff/ParallelRestrictedAdditiveSchwarz.hpp
  opm/autodiff/RateConverter.hpp
  opm/autodiff/RedistributeDataHandles.hpp
  opm/autodiff/SimFIBODetails.hpp
  opm/autodiff/SimulatorBase.hpp
  opm/autodiff/SimulatorBase_impl.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoil.hpp
  opm/autodiff/SimulatorSequentialBlackoil.hpp
  opm/autodiff/WellConnectionAuxiliaryModule.hpp
  opm/autodiff/WellDensitySegmented.hpp
  opm/autodiff/WellStateFullyImplicitBlackoil.hpp
  opm/autodiff/ThreadHandle.hpp
  opm/autodiff/VFPProperties.hpp
  opm/autodiff/VFPHelpers.hpp
  opm/autodiff/VFPProdProperties.hpp
  opm/autodiff/VFPInjProperties.hpp
  opm/autodiff/WellHelpers.hpp
  opm/autodiff/StandardWells.hpp
  opm/autodiff/StandardWells_impl.hpp
  opm/autodiff/WellInterface.hpp
  opm/autodiff/WellInterface_impl.hpp
  opm/autodiff/StandardWell.hpp
  opm/autodiff/StandardWell_impl.hpp
  opm/autodiff/MultisegmentWell.hpp
  opm/autodiff/MultisegmentWell_impl.hpp
  opm/autodiff/MSWellHelpers.hpp
  opm/autodiff/BlackoilWellModel.hpp
  opm/autodiff/BlackoilWellModel_impl.hpp
  opm/autodiff/MissingFeatures.hpp
  opm/core/linalg/LinearSolverPetsc.hpp
  opm/core/linalg/ParallelIstlInformation.hpp
  opm/core/pressure/legacy_well.h
  opm/core/pressure/tpfa/compr_quant_general.h
  opm/core/pressure/tpfa/compr_source.h
  opm/core/props/BlackoilPhases.hpp
  opm/core/props/BlackoilPropertiesInterface.hpp
  opm/core/props/IncompPropertiesInterface.hpp
  opm/core/props/IncompPropertiesShadow.hpp
  opm/core/props/IncompPropertiesShadow_impl.hpp
  opm/core/props/phaseUsageFromDeck.hpp
  opm/core/props/pvt/ThermalGasPvtWrapper.hpp
  opm/core/props/pvt/ThermalOilPvtWrapper.hpp
  opm/core/props/pvt/ThermalWaterPvtWrapper.hpp
  opm/core/props/rock/RockFromDeck.hpp
  opm/core/props/satfunc/RelpermDiagnostics.hpp
  opm/core/props/satfunc/SaturationPropsInterface.hpp
  opm/core/props/satfunc/RelpermDiagnostics_impl.hpp
  opm/core/simulator/BlackoilState.hpp
  opm/core/simulator/BlackoilStateToFluidState.hpp
  opm/core/simulator/EquilibrationHelpers.hpp
  opm/core/simulator/ExplicitArraysFluidState.hpp
  opm/core/simulator/ExplicitArraysSatDerivativesFluidState.hpp
  opm/core/simulator/SimulatorReport.hpp
  opm/core/simulator/WellState.hpp
  opm/core/simulator/initState.hpp
  opm/core/simulator/initStateEquil.hpp
  opm/core/simulator/initStateEquil_impl.hpp
  opm/core/simulator/initState_impl.hpp
  opm/core/utility/DataMap.hpp
  opm/core/utility/Event.hpp
  opm/core/utility/Event_impl.hpp
  opm/core/utility/initHydroCarbonState.hpp
  opm/core/utility/miscUtilities_impl.hpp
  opm/core/utility/share_obj.hpp
  opm/core/well_controls.h
  opm/core/wells.h
  opm/core/wells/InjectionSpecification.hpp
  opm/core/wells/ProductionSpecification.hpp
  opm/core/wells/WellCollection.hpp
  opm/core/wells/WellsGroup.hpp
  opm/core/wells/WellsManager.hpp
  opm/core/wells/DynamicListEconLimited.hpp
  opm/core/wells/WellsManager_impl.hpp
  opm/polymer/GravityColumnSolverPolymer.hpp
  opm/polymer/GravityColumnSolverPolymer_impl.hpp
  opm/polymer/IncompPropertiesDefaultPolymer.hpp
  opm/polymer/PolymerBlackoilState.hpp
  opm/polymer/SinglePointUpwindTwoPhasePolymer.hpp
  opm/polymer/Point2D.hpp
  opm/simulators/ParallelFileMerger.hpp
  opm/simulators/thresholdPressures.hpp
  opm/simulators/WellSwitchingLogger.hpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp
  opm/simulators/timestepping/AdaptiveTimeStepping.hpp
  opm/simulators/timestepping/AdaptiveTimeStepping_impl.hpp
  opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp
  opm/simulators/timestepping/ConvergenceReport.hpp
  opm/simulators/timestepping/TimeStepControl.hpp
  opm/simulators/timestepping/TimeStepControlInterface.hpp
  opm/simulators/timestepping/SimulatorTimer.hpp
  opm/simulators/timestepping/SimulatorTimerInterface.hpp
  opm/simulators/timestepping/gatherConvergenceReport.hpp
  )
