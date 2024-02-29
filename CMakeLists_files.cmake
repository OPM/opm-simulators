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
  ebos/eclmixingratecontrols.cc
  ebos/eclsolutioncontainers.cc
  ebos/ecltransmissibility.cc
  ebos/equil/equilibrationhelpers.cc
  ebos/equil/initstateequil.cc
  opm/core/props/BlackoilPhases.cpp
  opm/core/props/phaseUsageFromDeck.cpp
  opm/core/props/satfunc/RelpermDiagnostics.cpp
  opm/simulators/timestepping/SimulatorReport.cpp
  opm/simulators/flow/ActionHandler.cpp
  opm/simulators/flow/Banners.cpp
  opm/simulators/flow/ConvergenceOutputConfiguration.cpp
  opm/simulators/flow/ExtraConvergenceOutputThread.cpp
  opm/simulators/flow/FlowMain.cpp
  opm/simulators/flow/InterRegFlows.cpp
  opm/simulators/flow/KeywordValidation.cpp
  opm/simulators/flow/LogOutputHelper.cpp
  opm/simulators/flow/Main.cpp
  opm/simulators/flow/NonlinearSolver.cpp
  opm/simulators/flow/RSTConv.cpp
  opm/simulators/flow/RegionPhasePVAverage.cpp
  opm/simulators/flow/SimulatorReportBanners.cpp
  opm/simulators/flow/SimulatorSerializer.cpp
  opm/simulators/flow/ValidationFunctions.cpp
  opm/simulators/flow/partitionCells.cpp
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.cpp
  opm/simulators/linalg/FlexibleSolver1.cpp
  opm/simulators/linalg/FlexibleSolver2.cpp
  opm/simulators/linalg/FlexibleSolver3.cpp
  opm/simulators/linalg/FlexibleSolver4.cpp
  opm/simulators/linalg/FlexibleSolver5.cpp
  opm/simulators/linalg/FlexibleSolver6.cpp
  opm/simulators/linalg/ISTLSolver.cpp
  opm/simulators/linalg/MILU.cpp
  opm/simulators/linalg/ParallelIstlInformation.cpp
  opm/simulators/linalg/ParallelOverlappingILU0.cpp
  opm/simulators/linalg/PreconditionerFactory1.cpp
  opm/simulators/linalg/PreconditionerFactory2.cpp
  opm/simulators/linalg/PreconditionerFactory3.cpp
  opm/simulators/linalg/PreconditionerFactory4.cpp
  opm/simulators/linalg/PreconditionerFactory5.cpp
  opm/simulators/linalg/PreconditionerFactory6.cpp
  opm/simulators/linalg/PropertyTree.cpp
  opm/simulators/linalg/setupPropertyTree.cpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.cpp
  opm/simulators/timestepping/AdaptiveTimeStepping.cpp
  opm/simulators/timestepping/ConvergenceReport.cpp
  opm/simulators/timestepping/TimeStepControl.cpp
  opm/simulators/timestepping/SimulatorTimer.cpp
  opm/simulators/timestepping/SimulatorTimerInterface.cpp
  opm/simulators/timestepping/gatherConvergenceReport.cpp
  opm/simulators/utils/ComponentName.cpp
  opm/simulators/utils/compressPartition.cpp
  opm/simulators/utils/DeferredLogger.cpp
  opm/simulators/utils/gatherDeferredLogger.cpp
  opm/simulators/utils/ParallelFileMerger.cpp
  opm/simulators/utils/ParallelRestart.cpp
  opm/simulators/utils/PartiallySupportedFlowKeywords.cpp
  opm/simulators/utils/PressureAverage.cpp
  opm/simulators/utils/readDeck.cpp
  opm/simulators/utils/SerializationPackers.cpp
  opm/simulators/utils/UnsupportedFlowKeywords.cpp
  opm/simulators/wells/ALQState.cpp
  opm/simulators/wells/BlackoilWellModelConstraints.cpp
  opm/simulators/wells/BlackoilWellModelGeneric.cpp
  opm/simulators/wells/BlackoilWellModelGuideRates.cpp
  opm/simulators/wells/BlackoilWellModelRestart.cpp
  opm/simulators/wells/ConnFiltrateData.cpp
  opm/simulators/wells/FractionCalculator.cpp
  opm/simulators/wells/GasLiftCommon.cpp
  opm/simulators/wells/GasLiftGroupInfo.cpp
  opm/simulators/wells/GasLiftSingleWellGeneric.cpp
  opm/simulators/wells/GasLiftStage2.cpp
  opm/simulators/wells/GlobalWellInfo.cpp
  opm/simulators/wells/GroupEconomicLimitsChecker.cpp
  opm/simulators/wells/GroupState.cpp
  opm/simulators/wells/MSWellHelpers.cpp
  opm/simulators/wells/MultisegmentWellAssemble.cpp
  opm/simulators/wells/MultisegmentWellEquations.cpp
  opm/simulators/wells/MultisegmentWellEval.cpp
  opm/simulators/wells/MultisegmentWellGeneric.cpp
  opm/simulators/wells/MultisegmentWellPrimaryVariables.cpp
  opm/simulators/wells/MultisegmentWellSegments.cpp
  opm/simulators/wells/ParallelPAvgCalculator.cpp
  opm/simulators/wells/ParallelPAvgDynamicSourceData.cpp
  opm/simulators/wells/ParallelWBPCalculation.cpp
  opm/simulators/wells/ParallelWellInfo.cpp
  opm/simulators/wells/PerfData.cpp
  opm/simulators/wells/RateConverter.cpp
  opm/simulators/wells/SegmentState.cpp
  opm/simulators/wells/SingleWellState.cpp
  opm/simulators/wells/StandardWellAssemble.cpp
  opm/simulators/wells/StandardWellConnections.cpp
  opm/simulators/wells/StandardWellEquations.cpp
  opm/simulators/wells/StandardWellEval.cpp
  opm/simulators/wells/StandardWellPrimaryVariables.cpp
  opm/simulators/wells/TargetCalculator.cpp
  opm/simulators/wells/VFPHelpers.cpp
  opm/simulators/wells/VFPInjProperties.cpp
  opm/simulators/wells/VFPProdProperties.cpp
  opm/simulators/wells/WellAssemble.cpp
  opm/simulators/wells/WellBhpThpCalculator.cpp
  opm/simulators/wells/WellConnectionAuxiliaryModule.cpp
  opm/simulators/wells/WellConstraints.cpp
  opm/simulators/wells/WellConvergence.cpp
  opm/simulators/wells/WellFilterCake.cpp
  opm/simulators/wells/WellGroupConstraints.cpp
  opm/simulators/wells/WellGroupControls.cpp
  opm/simulators/wells/WellGroupHelpers.cpp
  opm/simulators/wells/WellHelpers.cpp
  opm/simulators/wells/WellInterfaceFluidSystem.cpp
  opm/simulators/wells/WellInterfaceGeneric.cpp
  opm/simulators/wells/WellInterfaceIndices.cpp
  opm/simulators/wells/WellProdIndexCalculator.cpp
  opm/simulators/wells/WellState.cpp
  opm/simulators/wells/WellTest.cpp
  opm/simulators/wells/WGState.cpp
  )


if (Damaris_FOUND AND MPI_FOUND AND USE_DAMARIS_LIB)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/utils/DamarisOutputModule.cpp
                                 opm/simulators/utils/DamarisKeywords.cpp
                                 opm/simulators/utils/initDamarisXmlFile.cpp
                                 ebos/damariswriter.cc
                                 opm/simulators/utils/DamarisVar.cpp
                                 opm/simulators/utils/GridDataOutput.cpp
  )
endif()
if(CUDA_FOUND)
  # CUISTL SOURCE
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/detail/CuBlasHandle.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.cu)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/detail/CuSparseHandle.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/CuVector.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/detail/vector_operations.cu)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/CuSparseMatrix.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/CuDILU.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/CuJac.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/CuSeqILU0.cpp)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/cuistl/set_device.cpp)

  # CUISTL HEADERS
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cuda_check_last_error.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/CuBlasHandle.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/CuSparseHandle.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/CuDILU.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/CuJac.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/CuVector.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/CuSparseMatrix.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/CuMatrixDescription.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/CuSparseResource.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/CuSparseResource_impl.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/safe_conversion.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cublas_wrapper.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cusparse_wrapper.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/cusparse_constants.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/vector_operations.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/has_function.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/preconditioner_should_call_post_pre.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/CuSeqILU0.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/detail/fix_zero_diagonal.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/PreconditionerConvertFieldTypeAdapter.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/CuOwnerOverlapCopy.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/SolverAdapter.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/CuBlockPreconditioner.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/PreconditionerHolder.hpp)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/linalg/cuistl/set_device.hpp)

endif()

if(USE_BDA_BRIDGE)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BdaBridge.cpp
                                 opm/simulators/linalg/bda/WellContributions.cpp
                                 opm/simulators/linalg/bda/MultisegmentWellContribution.cpp
                                 opm/simulators/linalg/ISTLSolverBda.cpp)
  if(OPENCL_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/BlockedMatrix.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/BILU0.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/Reorder.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/ChowPatelIlu.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/BISAI.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/CPR.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/opencl.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/openclKernels.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/OpenclMatrix.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/Preconditioner.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/openclSolverBackend.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/opencl/openclWellContributions.cpp)
  endif()
  if(ROCALUTION_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/rocalutionSolverBackend.cpp)
  endif()
  if(rocsparse_FOUND AND rocblas_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/rocsparseSolverBackend.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/rocsparseWellContributions.cpp)
  endif()
  if(CUDA_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/cuda/cusparseSolverBackend.cu)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/cuda/cuWellContributions.cu)
  endif()
  if(amgcl_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/amgclSolverBackend.cpp)
    if(CUDA_FOUND)
      list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/bda/cuda/amgclSolverBackend.cu)
    endif()
  endif()
endif()
if(MPI_FOUND)
  list(APPEND MAIN_SOURCE_FILES opm/simulators/utils/MPIPacker.cpp
                                opm/simulators/utils/ParallelEclipseState.cpp
                                opm/simulators/utils/ParallelNLDDPartitioningZoltan.cpp
                                opm/simulators/utils/ParallelSerialization.cpp
                                opm/simulators/utils/SetupZoltanParams.cpp)
  list(APPEND PUBLIC_HEADER_FILES opm/simulators/utils/MPIPacker.hpp
                                  opm/simulators/utils/MPISerializer.hpp)
endif()
if(HDF5_FOUND)
  list(APPEND MAIN_SOURCE_FILES opm/simulators/utils/HDF5File.cpp)
endif()

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/test_ALQState.cpp
  tests/test_aquifergridutils.cpp
  tests/test_blackoil_amg.cpp
  tests/test_convergenceoutputconfiguration.cpp
  tests/test_convergencereport.cpp
  tests/test_deferredlogger.cpp
  tests/test_dilu.cpp
  tests/test_equil.cc
  tests/test_extractMatrix.cpp
  tests/test_flexiblesolver.cpp
  tests/test_glift1.cpp
  tests/test_graphcoloring.cpp
  tests/test_GroupState.cpp
  tests/test_interregflows.cpp
  tests/test_invert.cpp
  tests/test_keyword_validator.cpp
  tests/test_LogOutputHelper.cpp
  tests/test_milu.cpp
  tests/test_multmatrixtransposed.cpp
  tests/test_norne_pvt.cpp
  tests/test_outputdir.cpp
  tests/test_parallel_wbp_sourcevalues.cpp
  tests/test_parallelwellinfo.cpp
  tests/test_partitionCells.cpp
  tests/test_preconditionerfactory.cpp
  tests/test_privarspacking.cpp
  tests/test_region_phase_pvaverage.cpp
  tests/test_relpermdiagnostics.cpp
  tests/test_RestartSerialization.cpp
  tests/test_rstconv.cpp
  tests/test_stoppedwells.cpp
  tests/test_timer.cpp
  tests/test_vfpproperties.cpp
  tests/test_wellmodel.cpp
  tests/test_wellprodindexcalculator.cpp
  tests/test_wellstate.cpp
  )

if (HAVE_ECL_INPUT)
  list(APPEND TEST_SOURCE_FILES tests/test_nonnc.cpp)
endif()

if(MPI_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_parallelistlinformation.cpp
                                tests/test_ParallelSerialization.cpp)
endif()

if(CUDA_FOUND)
  if(USE_BDA_BRIDGE)
    list(APPEND TEST_SOURCE_FILES tests/test_cusparseSolver.cpp)
  endif()
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_converttofloatadapter.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cublas_handle.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cublas_safe_call.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cusparse_safe_call.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cuda_safe_call.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cuda_check_last_error.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cujac.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cuowneroverlapcopy.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cuseqilu0.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cusparse_handle.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cuSparse_matrix_operations.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cusparsematrix.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cuvector.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_cuVector_operations.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_safe_conversion.cpp)
  list(APPEND TEST_SOURCE_FILES tests/cuistl/test_solver_adapter.cpp)


endif()

if(USE_BDA_BRIDGE)
  if(OPENCL_FOUND)
    list(APPEND TEST_SOURCE_FILES tests/test_openclSolver.cpp)
    list(APPEND TEST_SOURCE_FILES tests/test_solvetransposed3x3.cpp)
  list(APPEND TEST_SOURCE_FILES tests/test_csrToCscOffsetMap.cpp)
  endif()
  if(ROCALUTION_FOUND)
    list(APPEND TEST_SOURCE_FILES tests/test_rocalutionSolver.cpp)
  endif()
  if(rocsparse_FOUND AND rocblas_FOUND)
    list(APPEND TEST_SOURCE_FILES tests/test_rocsparseSolver.cpp)
  endif()
endif()

if(HDF5_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_HDF5File.cpp)
  list(APPEND TEST_SOURCE_FILES tests/test_HDF5Serializer.cpp)
endif()

list (APPEND TEST_DATA_FILES
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
  tests/equil_co2store_go.DATA
  tests/equil_co2store_gw.DATA
  tests/equil_wetgas.DATA
  tests/equil_liveoil.DATA
  tests/equil_humidwetgas.DATA
  tests/equil_rsvd_and_rvvd.DATA
  tests/equil_rsvd_and_rvvd_and_rvwvd.DATA
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
  tests/offset_map_matrix.txt
  tests/offset_map_matrix_transposed.txt
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
  tests/test10.partition
  )


# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  ebos/alucartesianindexmapper.hh
  ebos/collecttoiorank.hh
  ebos/collecttoiorank_impl.hh
  ebos/ebos.hh
  ebos/eclalugridvanguard.hh
  ebos/eclbaseaquifermodel.hh
  ebos/eclbasevanguard.hh
  ebos/eclcpgridvanguard.hh
  ebos/eclequilinitializer.hh
  ebos/eclfluxmodule.hh
  ebos/eclgenericcpgridvanguard.hh
  ebos/eclgenericoutputblackoilmodule.hh
  ebos/eclgenericproblem.hh
  ebos/eclgenericproblem_impl.hh
  ebos/eclgenericthresholdpressure.hh
  ebos/eclgenericthresholdpressure_impl.hh
  ebos/eclgenerictracermodel.hh
  ebos/eclgenerictracermodel_impl.hh
  ebos/eclgenericvanguard.hh
  ebos/eclgenericwriter.hh
  ebos/eclgenericwriter_impl.hh
  ebos/eclmixingratecontrols.hh
  ebos/eclnewtonmethod.hh
  ebos/ecloutputblackoilmodule.hh
  ebos/eclpolyhedralgridvanguard.hh
  ebos/eclproblem.hh
  ebos/eclproblem_properties.hh
  ebos/eclsolutioncontainers.hh
  ebos/eclthresholdpressure.hh
  ebos/ecltracermodel.hh
  ebos/ecltransmissibility.hh
  ebos/ecltransmissibility_impl.hh
  ebos/eclwriter.hh
  ebos/femcpgridcompat.hh
  ebos/vtkecltracermodule.hh
  opm/simulators/flow/ActionHandler.hpp
  opm/simulators/flow/Banners.hpp
  opm/simulators/flow/BlackoilModel.hpp
  opm/simulators/flow/BlackoilModelNldd.hpp
  opm/simulators/flow/BlackoilModelParameters.hpp
  opm/simulators/flow/ConvergenceOutputConfiguration.hpp
  opm/simulators/flow/countGlobalCells.hpp
  opm/simulators/flow/DummyGradientCalculator.hpp
  opm/simulators/flow/ExtraConvergenceOutputThread.hpp
  opm/simulators/flow/FlowMain.hpp
  opm/simulators/flow/FlowsData.hpp
  opm/simulators/flow/InterRegFlows.hpp
  opm/simulators/flow/KeywordValidation.hpp
  opm/simulators/flow/LogOutputHelper.hpp
  opm/simulators/flow/Main.hpp
  opm/simulators/flow/NonlinearSolver.hpp
  opm/simulators/flow/partitionCells.hpp
  opm/simulators/flow/priVarsPacking.hpp
  opm/simulators/flow/RSTConv.hpp
  opm/simulators/flow/RegionPhasePVAverage.hpp
  opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp
  opm/simulators/flow/SimulatorReportBanners.hpp
  opm/simulators/flow/SimulatorSerializer.hpp
  opm/simulators/flow/SubDomain.hpp
  opm/simulators/flow/ValidationFunctions.hpp
  opm/core/props/BlackoilPhases.hpp
  opm/core/props/phaseUsageFromDeck.hpp
  opm/core/props/satfunc/RelpermDiagnostics.hpp
  opm/simulators/timestepping/SimulatorReport.hpp
  opm/simulators/wells/SegmentState.hpp
  opm/simulators/wells/WellContainer.hpp
  opm/simulators/aquifers/AquiferAnalytical.hpp
  opm/simulators/aquifers/AquiferCarterTracy.hpp
  opm/simulators/aquifers/AquiferConstantFlux.hpp
  opm/simulators/aquifers/AquiferFetkovich.hpp
  opm/simulators/aquifers/AquiferGridUtils.hpp
  opm/simulators/aquifers/AquiferInterface.hpp
  opm/simulators/aquifers/AquiferNumerical.hpp
  opm/simulators/aquifers/BlackoilAquiferModel.hpp
  opm/simulators/aquifers/BlackoilAquiferModel_impl.hpp
  opm/simulators/linalg/bda/amgclSolverBackend.hpp
  opm/simulators/linalg/bda/BdaBridge.hpp
  opm/simulators/linalg/bda/BdaResult.hpp
  opm/simulators/linalg/bda/BdaSolver.hpp
  opm/simulators/linalg/bda/opencl/BILU0.hpp
  opm/simulators/linalg/bda/BlockedMatrix.hpp
  opm/simulators/linalg/bda/opencl/CPR.hpp
  opm/simulators/linalg/bda/cuda/cuda_header.hpp
  opm/simulators/linalg/bda/cuda/cusparseSolverBackend.hpp
  opm/simulators/linalg/bda/opencl/ChowPatelIlu.hpp
  opm/simulators/linalg/bda/opencl/BISAI.hpp
  opm/simulators/linalg/bda/Reorder.hpp
  opm/simulators/linalg/bda/opencl/opencl.hpp
  opm/simulators/linalg/bda/opencl/openclKernels.hpp
  opm/simulators/linalg/bda/opencl/OpenclMatrix.hpp
  opm/simulators/linalg/bda/opencl/Preconditioner.hpp
  opm/simulators/linalg/bda/opencl/openclSolverBackend.hpp
  opm/simulators/linalg/bda/opencl/openclWellContributions.hpp
  opm/simulators/linalg/bda/Matrix.hpp
  opm/simulators/linalg/bda/MultisegmentWellContribution.hpp
  opm/simulators/linalg/bda/rocalutionSolverBackend.hpp
  opm/simulators/linalg/bda/rocsparseSolverBackend.hpp
  opm/simulators/linalg/bda/rocsparseWellContributions.hpp
  opm/simulators/linalg/bda/WellContributions.hpp
  opm/simulators/linalg/amgcpr.hh
  opm/simulators/linalg/DILU.hpp
  opm/simulators/linalg/twolevelmethodcpr.hh
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp
  opm/simulators/linalg/ExtraSmoothers.hpp
  opm/simulators/linalg/FlexibleSolver.hpp
  opm/simulators/linalg/FlexibleSolver_impl.hpp
  opm/simulators/linalg/FlowLinearSolverParameters.hpp
  opm/simulators/linalg/GraphColoring.hpp
  opm/simulators/linalg/ISTLSolver.hpp
  opm/simulators/linalg/ISTLSolverBda.hpp
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
  opm/simulators/linalg/SmallDenseMatrixUtils.hpp
  opm/simulators/linalg/WellOperators.hpp
  opm/simulators/linalg/WriteSystemMatrixHelper.hpp
  opm/simulators/linalg/extractMatrix.hpp
  opm/simulators/linalg/findOverlapRowsAndColumns.hpp
  opm/simulators/linalg/getQuasiImpesWeights.hpp
  opm/simulators/linalg/setupPropertyTree.hpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp
  opm/simulators/timestepping/AdaptiveTimeStepping.hpp
  opm/simulators/timestepping/ConvergenceReport.hpp
  opm/simulators/timestepping/EclTimeSteppingParams.hpp
  opm/simulators/timestepping/TimeStepControl.hpp
  opm/simulators/timestepping/TimeStepControlInterface.hpp
  opm/simulators/timestepping/SimulatorTimer.hpp
  opm/simulators/timestepping/SimulatorTimerInterface.hpp
  opm/simulators/timestepping/gatherConvergenceReport.hpp
  opm/simulators/utils/ComponentName.hpp
  opm/simulators/utils/compressPartition.hpp
  opm/simulators/utils/ParallelFileMerger.hpp
  opm/simulators/utils/DeferredLoggingErrorHelpers.hpp
  opm/simulators/utils/DeferredLogger.hpp
  opm/simulators/utils/gatherDeferredLogger.hpp
  opm/simulators/utils/moduleVersion.hpp
  opm/simulators/utils/ParallelEclipseState.hpp
  opm/simulators/utils/ParallelNLDDPartitioningZoltan.hpp
  opm/simulators/utils/ParallelRestart.hpp
  opm/simulators/utils/PropsDataHandle.hpp
  opm/simulators/utils/SerializationPackers.hpp
  opm/simulators/utils/VectorVectorDataHandle.hpp
  opm/simulators/utils/PressureAverage.hpp
  opm/simulators/utils/readDeck.hpp
  opm/simulators/wells/ALQState.hpp
  opm/simulators/wells/BlackoilWellModel.hpp
  opm/simulators/wells/BlackoilWellModel_impl.hpp
  opm/simulators/wells/BlackoilWellModelConstraints.hpp
  opm/simulators/wells/BlackoilWellModelGeneric.hpp
  opm/simulators/wells/BlackoilWellModelGuideRates.hpp
  opm/simulators/wells/BlackoilWellModelRestart.hpp
  opm/simulators/wells/ConnFiltrateData.hpp
  opm/simulators/wells/FractionCalculator.hpp
  opm/simulators/wells/GasLiftCommon.hpp
  opm/simulators/wells/GasLiftGroupInfo.hpp
  opm/simulators/wells/GasLiftSingleWellGeneric.hpp
  opm/simulators/wells/GasLiftSingleWell.hpp
  opm/simulators/wells/GasLiftSingleWell_impl.hpp
  opm/simulators/wells/GasLiftStage2.hpp
  opm/simulators/wells/GasLiftWellState.hpp
  opm/simulators/wells/GlobalWellInfo.hpp
  opm/simulators/wells/GroupEconomicLimitsChecker.hpp
  opm/simulators/wells/GroupState.hpp
  opm/simulators/wells/MSWellHelpers.hpp
  opm/simulators/wells/MultisegmentWell.hpp
  opm/simulators/wells/MultisegmentWell_impl.hpp
  opm/simulators/wells/MultisegmentWellAssemble.hpp
  opm/simulators/wells/MultisegmentWellEquations.hpp
  opm/simulators/wells/MultisegmentWellEval.hpp
  opm/simulators/wells/MultisegmentWellGeneric.hpp
  opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp
  opm/simulators/wells/MultisegmentWellSegments.hpp
  opm/simulators/wells/ParallelPAvgCalculator.hpp
  opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp
  opm/simulators/wells/ParallelWBPCalculation.hpp
  opm/simulators/wells/ParallelWellInfo.hpp
  opm/simulators/wells/PerfData.hpp
  opm/simulators/wells/PerforationData.hpp
  opm/simulators/wells/RateConverter.hpp
  opm/simulators/wells/RegionAttributeHelpers.hpp
  opm/simulators/wells/RegionAverageCalculator.hpp
  opm/simulators/wells/SingleWellState.hpp
  opm/simulators/wells/StandardWell.hpp
  opm/simulators/wells/StandardWell_impl.hpp
  opm/simulators/wells/StandardWellAssemble.hpp
  opm/simulators/wells/StandardWellConnections.hpp
  opm/simulators/wells/StandardWellEquations.hpp
  opm/simulators/wells/StandardWellEval.hpp
  opm/simulators/wells/StandardWellPrimaryVariables.hpp
  opm/simulators/wells/TargetCalculator.hpp
  opm/simulators/wells/VFPHelpers.hpp
  opm/simulators/wells/VFPInjProperties.hpp
  opm/simulators/wells/VFPProdProperties.hpp
  opm/simulators/wells/VFPProperties.hpp
  opm/simulators/wells/WellAssemble.hpp
  opm/simulators/wells/WellBhpThpCalculator.hpp
  opm/simulators/wells/WellConnectionAuxiliaryModule.hpp
  opm/simulators/wells/WellConstraints.hpp
  opm/simulators/wells/WellConvergence.hpp
  opm/simulators/wells/WellFilterCake.hpp
  opm/simulators/wells/WellGroupConstraints.hpp
  opm/simulators/wells/WellGroupControls.hpp
  opm/simulators/wells/WellGroupHelpers.hpp
  opm/simulators/wells/WellHelpers.hpp
  opm/simulators/wells/WellInterface.hpp
  opm/simulators/wells/WellInterfaceGeneric.hpp
  opm/simulators/wells/WellInterface_impl.hpp
  opm/simulators/wells/WellProdIndexCalculator.hpp
  opm/simulators/wells/WellState.hpp
  opm/simulators/wells/WellTest.hpp
  opm/simulators/wells/WGState.hpp
  )

if (Damaris_FOUND AND MPI_FOUND AND USE_DAMARIS_LIB)
  list (APPEND PUBLIC_HEADER_FILES opm/simulators/utils/DamarisOutputModule.hpp
                                   opm/simulators/utils/DamarisKeywords.hpp
                                   ebos/damaris_properties.hh
                                   ebos/damariswriter.hh
                                   opm/simulators/utils/DamarisVar.hpp
                                   opm/simulators/utils/GridDataOutput.hpp
                                   opm/simulators/utils/GridDataOutput_impl.hpp
  )
endif()

if(HDF5_FOUND)
  list(APPEND PUBLIC_HEADER_FILES
    opm/simulators/utils/HDF5Serializer.hpp
    opm/simulators/utils/HDF5File.hpp
  )
  list(APPEND MAIN_SOURCE_FILES
    opm/simulators/utils/HDF5Serializer.cpp
  )
endif()

list (APPEND EXAMPLE_SOURCE_FILES
  examples/printvfp.cpp
)
if(HDF5_FOUND)
  list (APPEND EXAMPLE_SOURCE_FILES
    examples/opmrst_inspect.cpp
  )
endif()
