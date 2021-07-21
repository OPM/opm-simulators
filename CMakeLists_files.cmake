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
  ebos/collecttoiorank.cc
  ebos/eclgenericcpgridvanguard.cc
  ebos/eclgenericoutputblackoilmodule.cc
  ebos/eclgenericproblem.cc
  ebos/eclgenericthresholdpressure.cc
  ebos/eclgenerictracermodel.cc
  ebos/eclgenericvanguard.cc
  ebos/eclgenericwriter.cc
  ebos/ecltransmissibility.cc
  opm/core/props/phaseUsageFromDeck.cpp
  opm/core/props/satfunc/RelpermDiagnostics.cpp
  opm/simulators/timestepping/SimulatorReport.cpp
  opm/simulators/flow/countGlobalCells.cpp
  opm/simulators/flow/KeywordValidation.cpp
  opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.cpp
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.cpp
  opm/simulators/linalg/FlexibleSolver1.cpp
  opm/simulators/linalg/FlexibleSolver2.cpp
  opm/simulators/linalg/FlexibleSolver3.cpp
  opm/simulators/linalg/FlexibleSolver4.cpp
  opm/simulators/linalg/PropertyTree.cpp
  opm/simulators/linalg/setupPropertyTree.cpp
  opm/simulators/utils/PartiallySupportedFlowKeywords.cpp
  opm/simulators/utils/readDeck.cpp
  opm/simulators/utils/UnsupportedFlowKeywords.cpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.cpp
  opm/simulators/timestepping/AdaptiveTimeSteppingEbos.cpp
  opm/simulators/timestepping/TimeStepControl.cpp
  opm/simulators/timestepping/SimulatorTimer.cpp
  opm/simulators/timestepping/SimulatorTimerInterface.cpp
  opm/simulators/timestepping/gatherConvergenceReport.cpp
  opm/simulators/utils/DeferredLogger.cpp
  opm/simulators/utils/gatherDeferredLogger.cpp
  opm/simulators/utils/ParallelFileMerger.cpp
  opm/simulators/utils/ParallelRestart.cpp
  opm/simulators/wells/ALQState.cpp
  opm/simulators/wells/BlackoilWellModelGeneric.cpp
  opm/simulators/wells/GasLiftGroupInfo.cpp
  opm/simulators/wells/GasLiftSingleWellGeneric.cpp
  opm/simulators/wells/GasLiftStage2.cpp
  opm/simulators/wells/GlobalWellInfo.cpp
  opm/simulators/wells/GroupState.cpp
  opm/simulators/wells/MultisegmentWellEval.cpp
  opm/simulators/wells/MultisegmentWellGeneric.cpp
  opm/simulators/wells/ParallelWellInfo.cpp
  opm/simulators/wells/PerfData.cpp
  opm/simulators/wells/SegmentState.cpp
  opm/simulators/wells/StandardWellEval.cpp
  opm/simulators/wells/StandardWellGeneric.cpp
  opm/simulators/wells/TargetCalculator.cpp
  opm/simulators/wells/VFPHelpers.cpp
  opm/simulators/wells/VFPProdProperties.cpp
  opm/simulators/wells/VFPInjProperties.cpp
  opm/simulators/wells/WellGroupHelpers.cpp
  opm/simulators/wells/WellInterfaceEval.cpp
  opm/simulators/wells/WellInterfaceFluidSystem.cpp
  opm/simulators/wells/WellInterfaceGeneric.cpp
  opm/simulators/wells/WellInterfaceIndices.cpp
  opm/simulators/wells/WellProdIndexCalculator.cpp
  opm/simulators/wells/WellState.cpp
  opm/simulators/wells/WGState.cpp
  )

if(CUDA_FOUND)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/cusparseSolverBackend.cu)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/WellContributions.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/WellContributions.cu)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/MultisegmentWellContribution.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BdaBridge.cpp)
endif()
if(OPENCL_FOUND)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BlockedMatrix.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BILU0.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/Reorder.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/ChowPatelIlu.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/openclKernels.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/openclSolverBackend.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BdaBridge.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/WellContributions.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/MultisegmentWellContribution.cpp)
endif()
if(HAVE_FPGA)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BdaBridge.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/FPGAMatrix.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/FPGABILU0.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/FPGASolverBackend.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/FPGAUtils.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/MultisegmentWellContribution.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/WellContributions.cpp)
endif()
if(HAVE_AMGCL)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BdaBridge.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/WellContributions.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/MultisegmentWellContribution.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/amgclSolverBackend.cpp)
  if(CUDA_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/amgclSolverBackend.cu)
  endif()
endif()
if(MPI_FOUND)
  list(APPEND MAIN_SOURCE_FILES opm/simulators/utils/ParallelEclipseState.cpp
                                opm/simulators/utils/ParallelSerialization.cpp)
endif()

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/test_equil.cc
  tests/test_ecl_output.cc
  tests/test_blackoil_amg.cpp
  tests/test_convergencereport.cpp
  tests/test_flexiblesolver.cpp
  tests/test_preconditionerfactory.cpp
  tests/test_graphcoloring.cpp
  tests/test_vfpproperties.cpp
  tests/test_milu.cpp
  tests/test_multmatrixtransposed.cpp
  tests/test_wellmodel.cpp
  tests/test_deferredlogger.cpp
  tests/test_timer.cpp
  tests/test_invert.cpp
  tests/test_stoppedwells.cpp
  tests/test_relpermdiagnostics.cpp
  tests/test_norne_pvt.cpp
  tests/test_wellprodindexcalculator.cpp
  tests/test_wellstate.cpp
  tests/test_parallelwellinfo.cpp
  tests/test_glift1.cpp
  tests/test_keyword_validator.cpp
  tests/test_GroupState.cpp
  tests/test_ALQState.cpp
  )

if(MPI_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_parallelistlinformation.cpp
                                tests/test_ParallelRestart.cpp)
endif()
if(CUDA_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_cusparseSolver.cpp)
endif()
if(OPENCL_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_openclSolver.cpp)
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
  tests/matr33rep.txt
  tests/rhs3rep.txt
  tests/options_flexiblesolver.json
  tests/options_flexiblesolver_simple.json
  tests/GLIFT1.DATA
  tests/include/flowl_b_vfp.ecl
  tests/include/flowl_c_vfp.ecl
  tests/include/permx_model5.grdecl
  tests/include/pvt_live_oil_dgas.ecl
  tests/include/relperm.inc
  tests/include/rock.inc
  tests/include/summary.inc
  tests/include/test1_20x30x10.grdecl
  tests/include/well_vfp.ecl
  )


# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/simulators/flow/countGlobalCells.hpp
  opm/simulators/flow/BlackoilModelEbos.hpp
  opm/simulators/flow/BlackoilModelParametersEbos.hpp
  opm/simulators/flow/FlowMainEbos.hpp
  opm/simulators/flow/Main.hpp
  opm/simulators/flow/NonlinearSolverEbos.hpp
  opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp
  opm/simulators/flow/KeywordValidation.hpp
  opm/core/props/BlackoilPhases.hpp
  opm/core/props/phaseUsageFromDeck.hpp
  opm/core/props/satfunc/RelpermDiagnostics.hpp
  opm/simulators/timestepping/SimulatorReport.hpp
  opm/simulators/wells/SegmentState.hpp
  opm/simulators/wells/WellContainer.hpp
  opm/simulators/aquifers/AquiferInterface.hpp
  opm/simulators/aquifers/AquiferCarterTracy.hpp
  opm/simulators/aquifers/AquiferFetkovich.hpp
  opm/simulators/aquifers/AquiferNumerical.hpp
  opm/simulators/aquifers/BlackoilAquiferModel.hpp
  opm/simulators/aquifers/BlackoilAquiferModel_impl.hpp
  opm/simulators/linalg/bda/amgclSolverBackend.hpp
  opm/simulators/linalg/bda/BdaBridge.hpp
  opm/simulators/linalg/bda/BdaResult.hpp
  opm/simulators/linalg/bda/BdaSolver.hpp
  opm/simulators/linalg/bda/BILU0.hpp
  opm/simulators/linalg/bda/BlockedMatrix.hpp
  opm/simulators/linalg/bda/cuda_header.hpp
  opm/simulators/linalg/bda/cusparseSolverBackend.hpp
  opm/simulators/linalg/bda/ChowPatelIlu.hpp
  opm/simulators/linalg/bda/FPGAMatrix.hpp
  opm/simulators/linalg/bda/FPGABILU0.hpp
  opm/simulators/linalg/bda/FPGASolverBackend.hpp
  opm/simulators/linalg/bda/FPGAUtils.hpp
  opm/simulators/linalg/bda/Reorder.hpp
  opm/simulators/linalg/bda/ILUReorder.hpp
  opm/simulators/linalg/bda/opencl.hpp
  opm/simulators/linalg/bda/openclKernels.hpp
  opm/simulators/linalg/bda/openclSolverBackend.hpp
  opm/simulators/linalg/bda/MultisegmentWellContribution.hpp
  opm/simulators/linalg/bda/WellContributions.hpp
  opm/simulators/linalg/amgcpr.hh
  opm/simulators/linalg/twolevelmethodcpr.hh
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp
  opm/simulators/linalg/FlexibleSolver.hpp
  opm/simulators/linalg/FlexibleSolver_impl.hpp
  opm/simulators/linalg/FlowLinearSolverParameters.hpp
  opm/simulators/linalg/GraphColoring.hpp
  opm/simulators/linalg/ISTLSolverEbos.hpp
  opm/simulators/linalg/ISTLSolverEbosFlexible.hpp
  opm/simulators/linalg/MatrixBlock.hpp
  opm/simulators/linalg/MatrixMarketSpecializations.hpp
  opm/simulators/linalg/OwningBlockPreconditioner.hpp
  opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp
  opm/simulators/linalg/ParallelOverlappingILU0.hpp
  opm/simulators/linalg/ParallelRestrictedAdditiveSchwarz.hpp
  opm/simulators/linalg/ParallelIstlInformation.hpp
  opm/simulators/linalg/PressureSolverPolicy.hpp
  opm/simulators/linalg/PressureTransferPolicy.hpp
  opm/simulators/linalg/PreconditionerFactory.hpp
  opm/simulators/linalg/PreconditionerWithUpdate.hpp
  opm/simulators/linalg/PropertyTree.hpp
  opm/simulators/linalg/WellOperators.hpp
  opm/simulators/linalg/WriteSystemMatrixHelper.hpp
  opm/simulators/linalg/findOverlapRowsAndColumns.hpp
  opm/simulators/linalg/getQuasiImpesWeights.hpp
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
  opm/simulators/utils/ParallelEclipseState.hpp
  opm/simulators/utils/ParallelRestart.hpp
  opm/simulators/utils/PropsCentroidsDataHandle.hpp
  opm/simulators/wells/PerfData.hpp
  opm/simulators/wells/PerforationData.hpp
  opm/simulators/wells/RateConverter.hpp
  opm/simulators/utils/readDeck.hpp
  opm/simulators/wells/TargetCalculator.hpp
  opm/simulators/wells/WellConnectionAuxiliaryModule.hpp
  opm/simulators/wells/WellState.hpp
  opm/simulators/wells/GlobalWellInfo.hpp
  opm/simulators/wells/GroupState.hpp
  opm/simulators/wells/ALQState.hpp
  opm/simulators/wells/WGState.hpp
  opm/simulators/wells/VFPProperties.hpp
  opm/simulators/wells/VFPHelpers.hpp
  opm/simulators/wells/VFPInjProperties.hpp
  opm/simulators/wells/VFPProdProperties.hpp
  opm/simulators/wells/WellGroupHelpers.hpp
  opm/simulators/wells/WellHelpers.hpp
  opm/simulators/wells/WellInterface.hpp
  opm/simulators/wells/WellInterface_impl.hpp
  opm/simulators/wells/WellProdIndexCalculator.hpp
  opm/simulators/wells/StandardWell.hpp
  opm/simulators/wells/StandardWell_impl.hpp
  opm/simulators/wells/MultisegmentWell.hpp
  opm/simulators/wells/MultisegmentWell_impl.hpp
  opm/simulators/wells/MSWellHelpers.hpp
  opm/simulators/wells/BlackoilWellModel.hpp
  opm/simulators/wells/BlackoilWellModel_impl.hpp
  opm/simulators/wells/ParallelWellInfo.hpp
  )

list (APPEND EXAMPLE_SOURCE_FILES
  examples/printvfp.cpp
  )
