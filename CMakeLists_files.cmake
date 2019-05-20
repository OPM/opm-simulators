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
  ebos/nncsorter.cpp
  opm/autodiff/MPIUtilities.cpp
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
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.cpp
  opm/simulators/linalg/setupPropertyTree.cpp
  opm/simulators/timestepping/TimeStepControl.cpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.cpp
  opm/simulators/timestepping/SimulatorTimer.cpp
  opm/simulators/timestepping/gatherConvergenceReport.cpp
  opm/simulators/utils/DeferredLogger.cpp
  opm/simulators/utils/gatherDeferredLogger.cpp
  opm/simulators/utils/moduleVersion.cpp
  opm/simulators/wells/VFPProdProperties.cpp
  opm/simulators/wells/VFPInjProperties.cpp
  )

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/test_equil.cc
  tests/test_ecl_output.cc
  tests/test_blackoil_amg.cpp
  tests/test_convergencereport.cpp
  tests/test_flexiblesolver.cpp
  tests/test_graphcoloring.cpp
  tests/test_vfpproperties.cpp
  tests/test_milu.cpp
  tests/test_multmatrixtransposed.cpp
  tests/test_nncsorter.cpp
  tests/test_wellmodel.cpp
  tests/test_deferredlogger.cpp
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
  tests/test_wellstatefullyimplicitblackoil.cpp
  )

if(MPI_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_parallelistlinformation.cpp)
endif()

list (APPEND TEST_DATA_FILES
  tests/SUMMARY_DECK_NON_CONSTANT_POROSITY.DATA
  tests/equil_base.DATA
  tests/equil_capillary.DATA
  tests/equil_capillary_overlap.DATA
  tests/equil_capillary_swatinit.DATA
  tests/equil_deadfluids.DATA
  tests/equil_pbvd_and_pdvd.DATA
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
  tests/wells_no_perforation.data
  tests/matr33.txt
  tests/rhs3.txt
  tests/options_flexiblesolver.json
  )


# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/autodiff/BlackoilDetails.hpp
  opm/autodiff/BlackoilModelEbos.hpp
  opm/autodiff/BlackoilModelParametersEbos.hpp
  opm/autodiff/createGlobalCellArray.hpp
  opm/autodiff/FlowMainEbos.hpp
  opm/autodiff/GridInit.hpp
  opm/autodiff/IterationReport.hpp
  opm/autodiff/MPIUtilities.hpp
  opm/autodiff/NonlinearSolverEbos.hpp
  opm/autodiff/RateConverter.hpp
  opm/autodiff/SimFIBODetails.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilEbos.hpp
  opm/autodiff/MissingFeatures.hpp
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
  opm/core/wells/WellsManager_impl.hpp
  opm/simulators/aquifers/AquiferInterface.hpp
  opm/simulators/aquifers/AquiferCarterTracy.hpp
  opm/simulators/aquifers/AquiferFetkovich.hpp
  opm/simulators/aquifers/BlackoilAquiferModel.hpp
  opm/simulators/aquifers/BlackoilAquiferModel_impl.hpp
  opm/simulators/linalg/BlackoilAmg.hpp
  opm/simulators/linalg/BlackoilAmgCpr.hpp
  opm/simulators/linalg/amgcpr.hh
  opm/simulators/linalg/twolevelmethodcpr.hh
  opm/simulators/linalg/CPRPreconditioner.hpp
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp
  opm/simulators/linalg/FlexibleSolver.hpp
  opm/simulators/linalg/FlowLinearSolverParameters.hpp
  opm/simulators/linalg/GetQuasiImpesWeights.hpp
  opm/simulators/linalg/GraphColoring.hpp
  opm/simulators/linalg/ISTLSolverEbos.hpp
  opm/simulators/linalg/ISTLSolverEbosCpr.hpp
  opm/simulators/linalg/ISTLSolverEbosFlexible.hpp
  opm/simulators/linalg/MatrixBlock.hpp
  opm/simulators/linalg/MatrixMarketUtils.hpp
  opm/simulators/linalg/OwningBlockPreconditioner.hpp
  opm/simulators/linalg/OwningTwolevelPreconditioner.hpp
  opm/simulators/linalg/ParallelOverlappingILU0.hpp
  opm/simulators/linalg/ParallelRestrictedAdditiveSchwarz.hpp
  opm/simulators/linalg/ParallelIstlInformation.hpp
  opm/simulators/linalg/PressureSolverPolicy.hpp
  opm/simulators/linalg/PressureTransferPolicy.hpp
  opm/simulators/linalg/PreconditionerWithUpdate.hpp
  opm/simulators/linalg/makePreconditioner.hpp
  opm/simulators/linalg/setupPropertyTree.hpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp
  opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp
  opm/simulators/timestepping/ConvergenceReport.hpp
  opm/simulators/timestepping/TimeStepControl.hpp
  opm/simulators/timestepping/TimeStepControlInterface.hpp
  opm/simulators/timestepping/SimulatorTimer.hpp
  opm/simulators/timestepping/SimulatorTimerInterface.hpp
  opm/simulators/timestepping/gatherConvergenceReport.hpp
  opm/simulators/utils/ParallelFileMerger.hpp
  opm/simulators/utils/DeferredLoggingErrorHelpers.hpp
  opm/simulators/utils/DeferredLogger.hpp
  opm/simulators/utils/gatherDeferredLogger.hpp
  opm/simulators/utils/moduleVersion.hpp
  opm/simulators/wells/WellConnectionAuxiliaryModule.hpp
  opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp
  opm/simulators/wells/VFPProperties.hpp
  opm/simulators/wells/VFPHelpers.hpp
  opm/simulators/wells/VFPInjProperties.hpp
  opm/simulators/wells/VFPProdProperties.hpp
  opm/simulators/wells/WellHelpers.hpp
  opm/simulators/wells/WellInterface.hpp
  opm/simulators/wells/WellInterface_impl.hpp
  opm/simulators/wells/StandardWell.hpp
  opm/simulators/wells/StandardWell_impl.hpp
  opm/simulators/wells/MultisegmentWell.hpp
  opm/simulators/wells/MultisegmentWell_impl.hpp
  opm/simulators/wells/StandardWellV.hpp
  opm/simulators/wells/StandardWellV_impl.hpp
  opm/simulators/wells/MSWellHelpers.hpp
  opm/simulators/wells/BlackoilWellModel.hpp
  opm/simulators/wells/BlackoilWellModel_impl.hpp
  )
