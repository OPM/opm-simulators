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
  opm/autodiff/ExtractParallelGridInformationToISTL.cpp
  opm/autodiff/moduleVersion.cpp
  opm/autodiff/MPIUtilities.cpp
  opm/autodiff/VFPProdProperties.cpp
  opm/autodiff/VFPInjProperties.cpp
  opm/autodiff/MissingFeatures.cpp
  opm/core/props/rock/RockFromDeck.cpp
  opm/core/props/satfunc/RelpermDiagnostics.cpp
  opm/core/simulator/SimulatorReport.cpp
  opm/core/wells/InjectionSpecification.cpp
  opm/core/wells/ProductionSpecification.cpp
  opm/core/wells/WellCollection.cpp
  opm/core/wells/WellsGroup.cpp
  opm/core/wells/WellsManager.cpp
  opm/core/wells/well_controls.c
  opm/core/wells/wells.c
  opm/simulators/WellSwitchingLogger.cpp
  opm/simulators/timestepping/TimeStepControl.cpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.cpp
  opm/simulators/timestepping/SimulatorTimer.cpp
  opm/simulators/timestepping/gatherConvergenceReport.cpp
  )

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/test_blackoil_amg.cpp
  tests/test_convergencereport.cpp
  tests/test_graphcoloring.cpp
  tests/test_vfpproperties.cpp
  tests/test_milu.cpp
  tests/test_multmatrixtransposed.cpp
  tests/test_wellmodel.cpp
  tests/test_wellswitchlogger.cpp
  tests/test_timer.cpp
  tests/test_invert.cpp
  tests/test_wells.cpp
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
  tests/wells_manager_data.data
  tests/wells_manager_data_expanded.data
  tests/wells_manager_data_wellSTOP.data
  tests/wells_group.data
  tests/wells_stopped.data
  tests/relpermDiagnostics.DATA
  tests/norne_pvt.data
  )


# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/autodiff/AquiferCarterTracy.hpp
  opm/autodiff/BlackoilAmg.hpp
  opm/autodiff/BlackoilDetails.hpp
  opm/autodiff/BlackoilModelParametersEbos.hpp
  opm/autodiff/BlackoilAquiferModel.hpp
  opm/autodiff/BlackoilAquiferModel_impl.hpp
  opm/autodiff/CPRPreconditioner.hpp
  opm/autodiff/createGlobalCellArray.hpp
  opm/autodiff/ExtractParallelGridInformationToISTL.hpp
  opm/autodiff/FlowLinearSolverParameters.hpp
  opm/autodiff/FlowMainEbos.hpp
  opm/autodiff/GraphColoring.hpp
  opm/autodiff/ISTLSolverEbos.hpp
  opm/autodiff/IterationReport.hpp
  opm/autodiff/MatrixBlock.hpp
  opm/autodiff/moduleVersion.hpp
  opm/autodiff/MPIUtilities.hpp
  opm/autodiff/NonlinearSolverEbos.hpp
  opm/autodiff/ParallelOverlappingILU0.hpp
  opm/autodiff/ParallelRestrictedAdditiveSchwarz.hpp
  opm/autodiff/RateConverter.hpp
  opm/autodiff/SimFIBODetails.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp
  opm/autodiff/WellConnectionAuxiliaryModule.hpp
  opm/autodiff/WellStateFullyImplicitBlackoil.hpp
  opm/autodiff/VFPProperties.hpp
  opm/autodiff/VFPHelpers.hpp
  opm/autodiff/VFPInjProperties.hpp
  opm/autodiff/VFPProdProperties.hpp
  opm/autodiff/WellHelpers.hpp
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
  opm/core/linalg/ParallelIstlInformation.hpp
  opm/core/props/BlackoilPhases.hpp
  opm/core/props/phaseUsageFromDeck.hpp
  opm/core/props/rock/RockFromDeck.hpp
  opm/core/props/satfunc/RelpermDiagnostics.hpp
  opm/core/props/satfunc/RelpermDiagnostics_impl.hpp
  opm/core/simulator/SimulatorReport.hpp
  opm/core/simulator/WellState.hpp
  opm/core/well_controls.h
  opm/core/wells.h
  opm/core/wells/InjectionSpecification.hpp
  opm/core/wells/ProductionSpecification.hpp
  opm/core/wells/WellCollection.hpp
  opm/core/wells/WellsGroup.hpp
  opm/core/wells/WellsManager.hpp
  opm/core/wells/DynamicListEconLimited.hpp
  opm/core/wells/WellsManager_impl.hpp
  opm/simulators/ParallelFileMerger.hpp
  opm/simulators/WellSwitchingLogger.hpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp
  opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp
  opm/simulators/timestepping/ConvergenceReport.hpp
  opm/simulators/timestepping/TimeStepControl.hpp
  opm/simulators/timestepping/TimeStepControlInterface.hpp
  opm/simulators/timestepping/SimulatorTimer.hpp
  opm/simulators/timestepping/SimulatorTimerInterface.hpp
  opm/simulators/timestepping/gatherConvergenceReport.hpp
  )
