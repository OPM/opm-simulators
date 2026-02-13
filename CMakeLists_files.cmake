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

# This macro adds a cuda/hip source file to the correct source file list
# it takes in the list to add it to, the path to the cuistl directory, and then
# the rest of the file path after cuistl. The reason for splitting this into to
# paths is to simplify replacing the cuistl part with hipistl.
# Cuda files are added as they are, whereas hip files should be added after
# hipification, we a dependency that will trigger when the cuda source code is
# changed.
macro (ADD_CUDA_OR_HIP_FILE LIST DIR FILE)
  set (cuda_file_path "${PROJECT_SOURCE_DIR}/${DIR}/gpuistl/${FILE}")

  if(CUDA_FOUND)
    list (APPEND ${LIST} "${DIR}/gpuistl/${FILE}")
  elseif(CONVERT_CUDA_TO_HIP)
    # we must hipify the code
    # and include the correct path which is in the build/binary dir
    string(REPLACE ".cu" ".hip" HIP_SOURCE_FILE ${FILE})
    set (hip_file_path "${PROJECT_BINARY_DIR}/${DIR}/gpuistl_hip/${HIP_SOURCE_FILE}")

    # add a custom command that will hipify
    add_custom_command(
        OUTPUT
          ${hip_file_path}
        COMMAND
          bash
          ${PROJECT_SOURCE_DIR}/bin/hipify_file.sh
          ${cuda_file_path}
          ${hip_file_path}
          $<TARGET_FILE:hipify-perl>
        DEPENDS
          ${cuda_file_path}
        COMMENT
          "Rehipifying because of change in ${cuda_file_path}"
    )

    # set_source_files_properties(${relpath} PROPERTIES LANGUAGE HIP)
    if("${LIST}" STREQUAL "PUBLIC_HEADER_FILES")
      file(RELATIVE_PATH relpath ${PROJECT_BINARY_DIR} ${hip_file_path})
      list(APPEND GENERATED_HEADER_FILES ${relpath})
    else()
      file(RELATIVE_PATH relpath ${PROJECT_SOURCE_DIR} ${hip_file_path})
      list(APPEND ${LIST} ${relpath})
    endif()
    list(APPEND ${LIST}_HIPIFIED ${hip_file_path})
  endif()
endmacro()

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
  flowexperimental/BlackOilEnergyIntensiveQuantitiesGlobalIndex.hpp
  flowexperimental/BlackOilIntensiveQuantitiesGlobalIndex.hpp
  flowexperimental/comp/EmptyModel.hpp
  flowexperimental/comp/flowexp_comp.hpp
  flowexperimental/comp/wells/CompWellModel.hpp
  flowexperimental/comp/wells/CompWellModel_impl.hpp
  flowexperimental/comp/wells/CompWellEquations.hpp
  flowexperimental/comp/wells/CompWellEquations_impl.hpp
  flowexperimental/comp/wells/CompWell.hpp
  flowexperimental/comp/wells/CompWell_impl.hpp
  flowexperimental/comp/wells/CompWellInterface.hpp
  flowexperimental/comp/wells/CompWellInterface_impl.hpp
  flowexperimental/comp/wells/CompWellPrimaryVariables.hpp
  flowexperimental/comp/wells/CompWellPrimaryVariables_impl.hpp
  flowexperimental/comp/wells/CompWellState.hpp
  flowexperimental/comp/wells/CompWellState_impl.hpp
  flowexperimental/comp/wells/SingleCompWellState.hpp
  flowexperimental/comp/wells/SingleCompWellState_impl.hpp
  flowexperimental/FIBlackOilModelNoCache.hpp
  flowexperimental/flowexp.hpp
  flowexperimental/FlowExpNewtonMethod.hpp
  opm/models/blackoil/blackoilbioeffectsparams.cpp
  opm/models/blackoil/blackoilbrineparams.cpp
  opm/models/blackoil/blackoilextboparams.cpp
  opm/models/blackoil/blackoilfoamparams.cpp
  opm/models/blackoil/blackoilnewtonmethodparams.cpp
  opm/models/blackoil/blackoilpolymerparams.cpp
  opm/models/blackoil/blackoilsolventparams.cpp
  opm/models/io/vtkblackoilbioeffectsparams.cpp
  opm/models/io/vtkblackoilenergyparams.cpp
  opm/models/io/vtkblackoilpolymerparams.cpp
  opm/models/io/vtkblackoilparams.cpp
  opm/models/io/vtkblackoilsolventparams.cpp
  opm/models/io/vtkcompositionparams.cpp
  opm/models/io/vtkdiffusionparams.cpp
  opm/models/io/vtkdiscretefractureparams.cpp
  opm/models/io/vtkenergyparams.cpp
  opm/models/io/vtkmultiphaseparams.cpp
  opm/models/io/vtkphasepresenceparams.cpp
  opm/models/io/vtkprimaryvarsparams.cpp
  opm/models/io/vtkptflashparams.cpp
  opm/models/io/vtktemperatureparams.cpp
  opm/models/io/restart.cpp
  opm/models/nonlinear/newtonmethodparams.cpp
  opm/models/parallel/tasklets.cpp
  opm/models/parallel/threadmanager.cpp
  opm/models/utils/parametersystem.cpp
  opm/models/utils/simulatorutils.cpp
  opm/models/utils/terminal.cpp
  opm/models/utils/timer.cpp
  opm/simulators/flow/ActionHandler.cpp
  opm/simulators/flow/Banners.cpp
  opm/simulators/flow/BioeffectsContainer.cpp
  opm/simulators/flow/BlackoilModelParameters.cpp
  opm/simulators/flow/BlackoilModelConvergenceMonitor.cpp
  opm/simulators/flow/CO2H2Container.cpp
  opm/simulators/flow/CollectDataOnIORank.cpp
  opm/simulators/flow/CompositionalContainer.cpp
  opm/simulators/flow/ConvergenceOutputConfiguration.cpp
  opm/simulators/flow/EclGenericWriter.cpp
  opm/simulators/flow/ExtboContainer.cpp
  opm/simulators/flow/ExtraConvergenceOutputThread.cpp
  opm/simulators/flow/FIPContainer.cpp
  opm/simulators/flow/FlowGenericProblem.cpp
  opm/simulators/flow/FlowGenericVanguard.cpp
  opm/simulators/flow/FlowProblemParameters.cpp
  opm/simulators/flow/FlowsContainer.cpp
  opm/simulators/flow/FlowUtils.cpp
  opm/simulators/flow/GenericCpGridVanguard.cpp
  opm/simulators/flow/GenericOutputBlackoilModule.cpp
  opm/simulators/flow/GenericTemperatureModel.cpp
  opm/simulators/flow/GenericThresholdPressure.cpp
  opm/simulators/flow/GenericTracerModel.cpp
  opm/simulators/flow/HybridNewtonConfig.cpp
  opm/simulators/flow/InterRegFlows.cpp
  opm/simulators/flow/KeywordValidation.cpp
  opm/simulators/flow/LogOutputHelper.cpp
  opm/simulators/flow/Main.cpp
  opm/simulators/flow/MechContainer.cpp
  opm/simulators/flow/MixingRateControls.cpp
  opm/simulators/flow/NlddReporting.cpp
  opm/simulators/flow/NonlinearSolver.cpp
  opm/simulators/flow/partitionCells.cpp
  opm/simulators/flow/RFTContainer.cpp
  opm/simulators/flow/RSTConv.cpp
  opm/simulators/flow/RegionPhasePVAverage.cpp
  opm/simulators/flow/SimulatorConvergenceOutput.cpp
  opm/simulators/flow/SimulatorFullyImplicitBlackoil.cpp
  opm/simulators/flow/SimulatorReportBanners.cpp
  opm/simulators/flow/SimulatorSerializer.cpp
  opm/simulators/flow/SolutionContainers.cpp
  opm/simulators/flow/TracerContainer.cpp
  opm/simulators/flow/Transmissibility.cpp
  opm/simulators/flow/ValidationFunctions.cpp
  opm/simulators/flow/equil/EquilibrationHelpers.cpp
  opm/simulators/flow/equil/InitStateEquil.cpp
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.cpp
  opm/simulators/linalg/FlexibleSolver1.cpp
  opm/simulators/linalg/FlexibleSolver2.cpp
  opm/simulators/linalg/FlexibleSolver3.cpp
  opm/simulators/linalg/FlexibleSolver4.cpp
  opm/simulators/linalg/FlexibleSolver5.cpp
  opm/simulators/linalg/FlexibleSolver6.cpp
  opm/simulators/linalg/FlexibleSolver7.cpp
  opm/simulators/linalg/FlowLinearSolverParameters.cpp
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
  opm/simulators/linalg/PreconditionerFactory7.cpp
  opm/simulators/linalg/PropertyTree.cpp
  opm/simulators/linalg/setupPropertyTree.cpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.cpp
  opm/simulators/timestepping/AdaptiveTimeStepping.cpp
  opm/simulators/timestepping/ConvergenceReport.cpp
  opm/simulators/timestepping/EclTimeSteppingParams.cpp
  opm/simulators/timestepping/SimulatorReport.cpp
  opm/simulators/timestepping/SimulatorTimer.cpp
  opm/simulators/timestepping/SimulatorTimerInterface.cpp
  opm/simulators/timestepping/TimeStepControl.cpp
  opm/simulators/timestepping/gatherConvergenceReport.cpp
  opm/simulators/utils/ComponentName.cpp
  opm/simulators/utils/DeferredLogger.cpp
  opm/simulators/utils/FullySupportedFlowKeywords.cpp
  opm/simulators/utils/InstantiationIndicesMacros.hpp
  opm/simulators/utils/ParallelFileMerger.cpp
  opm/simulators/utils/ParallelRestart.cpp
  opm/simulators/utils/PartiallySupportedFlowKeywords.cpp
  opm/simulators/utils/PressureAverage.cpp
  opm/simulators/utils/SerializationPackers.cpp
  opm/simulators/utils/UnsupportedFlowKeywords.cpp
  opm/simulators/utils/compressPartition.cpp
  opm/simulators/utils/gatherDeferredLogger.cpp
  opm/simulators/utils/readDeck.cpp
  opm/simulators/utils/satfunc/RelpermDiagnostics.cpp
  opm/simulators/wells/ALQState.cpp
  opm/simulators/wells/BlackoilWellModelConstraints.cpp
  opm/simulators/wells/BlackoilWellModelGasLift.cpp
  opm/simulators/wells/BlackoilWellModelGeneric.cpp
  opm/simulators/wells/BlackoilWellModelGuideRates.cpp
  opm/simulators/wells/BlackoilWellModelNetworkGeneric.cpp
  opm/simulators/wells/BlackoilWellModelNldd.cpp
  opm/simulators/wells/BlackoilWellModelRestart.cpp
  opm/simulators/wells/BlackoilWellModelWBP.cpp
  opm/simulators/wells/ConnFiltrateData.cpp
  opm/simulators/wells/GuideRateHandler.cpp
  opm/simulators/wells/FractionCalculator.cpp
  opm/simulators/wells/GasLiftCommon.cpp
  opm/simulators/wells/GasLiftGroupInfo.cpp
  opm/simulators/wells/GasLiftSingleWellGeneric.cpp
  opm/simulators/wells/GasLiftStage2.cpp
  opm/simulators/wells/GlobalWellInfo.cpp
  opm/simulators/wells/GroupEconomicLimitsChecker.cpp
  opm/simulators/wells/GroupState.cpp
  opm/simulators/wells/GroupStateHelper.cpp
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
  opm/simulators/wells/RatioCalculator.cpp
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
  opm/simulators/wells/WellConstraints.cpp
  opm/simulators/wells/WellConvergence.cpp
  opm/simulators/wells/WellFilterCake.cpp
  opm/simulators/wells/WellGroupConstraints.cpp
  opm/simulators/wells/WellGroupControls.cpp
  opm/simulators/wells/WellHelpers.cpp
  opm/simulators/wells/WellInterfaceFluidSystem.cpp
  opm/simulators/wells/WellInterfaceGeneric.cpp
  opm/simulators/wells/WellInterfaceIndices.cpp
  opm/simulators/wells/WellProdIndexCalculator.cpp
  opm/simulators/wells/WellState.cpp
  opm/simulators/wells/WellTest.cpp
  opm/simulators/wells/WGState.cpp
  )

if (HAVE_AVX2_EXTENSION)
  set (AVX2_SOURCE_FILES
    opm/simulators/linalg/mixed/bsr.c
    opm/simulators/linalg/mixed/prec.c
    opm/simulators/linalg/mixed/bslv.c)
  list (APPEND MAIN_SOURCE_FILES
    ${AVX2_SOURCE_FILES})
endif()

if (HAVE_ECL_INPUT)
  list (APPEND MAIN_SOURCE_FILES
    opm/simulators/utils/satfunc/GasPhaseConsistencyChecks.cpp
    opm/simulators/utils/satfunc/OilPhaseConsistencyChecks.cpp
    opm/simulators/utils/satfunc/PhaseCheckBase.cpp
    opm/simulators/utils/satfunc/SatfuncConsistencyCheckManager.cpp
    opm/simulators/utils/satfunc/SatfuncConsistencyChecks.cpp
    opm/simulators/utils/satfunc/ScaledSatfuncCheckPoint.cpp
    opm/simulators/utils/satfunc/ThreePointHorizontalConsistencyChecks.cpp
    opm/simulators/utils/satfunc/UnscaledSatfuncCheckPoint.cpp
    opm/simulators/utils/satfunc/WaterPhaseConsistencyChecks.cpp
  )
endif()

if (Damaris_FOUND AND MPI_FOUND AND USE_DAMARIS_LIB)
  list (APPEND MAIN_SOURCE_FILES
    opm/simulators/flow/DamarisParameters.cpp
    opm/simulators/flow/DamarisWriter.cpp
    opm/simulators/utils/DamarisKeywords.cpp
    opm/simulators/utils/DamarisOutputModule.cpp
    opm/simulators/utils/DamarisVar.cpp
    opm/simulators/utils/GridDataOutput.cpp
    opm/simulators/utils/initDamarisXmlFile.cpp
  )
endif()

# add these files if we should compile the hip code
if(CUDA_FOUND OR hip_FOUND)
  list(APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpuistl/device_management.hpp) # should not be hipified to make main independant of library
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg device_management.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/CuBlasHandle.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/gpusparse_matrix_operations.cu)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/CuSparseHandle.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/preconditionerKernels/DILUKernels.cu)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/preconditionerKernels/ILU_variants_helper_kernels.cu)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/preconditionerKernels/ILU0Kernels.cu)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/preconditionerKernels/JacKernels.cu)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/kernel_enums.hpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg GpuVector.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg GpuView.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/cpr_amg_operations.cu)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/vector_operations.cu)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg GpuSparseMatrix.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg GpuSparseMatrixGeneric.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg GpuDILU.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg OpmGpuILU0.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg GpuJac.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg GpuSeqILU0.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg set_device.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg detail/FlexibleSolverWrapper.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg FlexibleSolver_gpu_instantiate.cpp)
  ADD_CUDA_OR_HIP_FILE(MAIN_SOURCE_FILES opm/simulators/linalg PreconditionerFactory_gpu_instantiate.cpp)


  # HEADERS
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/autotuner.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/coloringAndReorderingUtils.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/gpu_safe_call.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/gpusparse_matrix_operations.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/gpusparse_matrix_utilities.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/cusparse_safe_call.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/cublas_safe_call.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/cuda_check_last_error.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/CuBlasHandle.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/CuSparseHandle.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuBuffer.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/preconditionerKernels/DILUKernels.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/preconditionerKernels/ILU_variants_helper_kernels.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/preconditionerKernels/ILU0Kernels.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/preconditionerKernels/JacKernels.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuDILU.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg OpmGpuILU0.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuJac.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuVector.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuView.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuSparseMatrix.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuSparseMatrixWrapper.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuSparseMatrixGeneric.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg MiniMatrix.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg MiniVector.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/CuMatrixDescription.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/CuSparseResource.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/CuSparseResource_impl.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/safe_conversion.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/gpu_type_detection.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/cublas_wrapper.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/cusparse_wrapper.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/gpu_constants.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/vector_operations.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/cpr_amg_operations.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/has_function.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/preconditioner_should_call_post_pre.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/deviceBlockOperations.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/gpuThreadUtils.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg PreconditionerAdapter.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuSeqILU0.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/fix_zero_diagonal.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg PreconditionerConvertFieldTypeAdapter.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg SolverAdapter.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuBlockPreconditioner.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg PreconditionerHolder.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg set_device.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg gpu_smart_pointer.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg gpu_resources.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/is_gpu_pointer.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg PreconditionerCPUMatrixToGPUMatrix.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg ISTLSolverGPUISTL.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/FlexibleSolverWrapper.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg AmgxInterface.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg HypreInterface.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg hypreinterface/HypreCpuTransfers.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg hypreinterface/HypreDataStructures.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg hypreinterface/HypreErrorHandling.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg hypreinterface/HypreGpuTransfers.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg hypreinterface/HypreSetup.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg hypreinterface/HypreUtils.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg PinnedMemoryHolder.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuPressureTransferPolicy.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg detail/gpu_preconditioner_utils.hpp)
  ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg MiniVector.hpp)

  if(MPI_FOUND)
    ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuOwnerOverlapCopy.hpp)
    ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuSender.hpp)
    ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuObliviousMPISender.hpp)
    ADD_CUDA_OR_HIP_FILE(PUBLIC_HEADER_FILES opm/simulators/linalg GpuAwareMPISender.hpp)
  endif()
endif()

if(USE_GPU_BRIDGE)
  list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/GpuBridge.cpp
                                 opm/simulators/linalg/gpubridge/CprCreation.cpp
                                 opm/simulators/linalg/gpubridge/Misc.cpp
                                 opm/simulators/linalg/gpubridge/WellContributions.cpp
                                 opm/simulators/linalg/gpubridge/MultisegmentWellContribution.cpp
                                 opm/simulators/linalg/ISTLSolverGpuBridge.cpp)
  if(OPENCL_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/BlockedMatrix.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/openclBILU0.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/Reorder.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/ChowPatelIlu.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/openclBISAI.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/openclCPR.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/opencl.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/openclKernels.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/OpenclMatrix.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/openclPreconditioner.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/openclSolverBackend.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/opencl/openclWellContributions.cpp)
  endif()
  if(ROCALUTION_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/rocalutionSolverBackend.cpp)
  endif()
  if(rocsparse_FOUND AND rocblas_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/rocsparseCPR.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/rocsparseBILU0.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/rocsparsePreconditioner.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/rocsparseSolverBackend.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/rocsparseWellContributions.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/hipKernels.cpp)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/rocm/rocsparseMatrix.cpp)
  endif()
  if(CUDA_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/cuda/cusparseSolverBackend.cu)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/cuda/cuWellContributions.cu)
  endif()
  if(amgcl_FOUND)
    list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/amgclSolverBackend.cpp)
    if(CUDA_FOUND)
      list (APPEND MAIN_SOURCE_FILES opm/simulators/linalg/gpubridge/cuda/amgclSolverBackend.cu)
    endif()
  endif()
endif()
if(MPI_FOUND)
  list(APPEND MAIN_SOURCE_FILES opm/simulators/utils/MPIPacker.cpp
                                opm/simulators/utils/ParallelEclipseState.cpp
                                opm/simulators/utils/ParallelNLDDPartitioningZoltan.cpp
                                opm/simulators/utils/ParallelSerialization.cpp
                                opm/simulators/utils/SetupPartitioningParams.cpp)
  list(APPEND PUBLIC_HEADER_FILES opm/simulators/utils/MPIPacker.hpp
                                  opm/simulators/utils/MPISerializer.hpp)
endif()
if(HDF5_FOUND)
  list(APPEND MAIN_SOURCE_FILES opm/simulators/utils/HDF5File.cpp)
endif()

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/models/test_quadrature.cpp
  tests/models/test_propertysystem.cpp
  tests/models/test_tasklets.cpp
  tests/models/test_tasklets_failure.cpp
  tests/test_ALQState.cpp
  tests/test_aquifergridutils.cpp
  tests/test_blackoil_amg.cpp
  tests/test_convergenceoutputconfiguration.cpp
  tests/test_convergencereport.cpp
  tests/test_deferredlogger.cpp
  tests/test_dilu.cpp
  tests/test_group_higher_constraints.cpp
  tests/test_equil.cpp
  tests/test_extractMatrix.cpp
  tests/test_flexiblesolver.cpp
  tests/test_gconsump.cpp
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
  tests/test_parametersystem.cpp
  tests/test_parallel_wbp_sourcevalues.cpp
  tests/test_parallelwellinfo.cpp
  tests/test_partitionCells.cpp
  tests/test_preconditionerfactory.cpp
  tests/test_privarspacking.cpp
  tests/test_propertytree.cpp
  tests/test_region_phase_pvaverage.cpp
  tests/test_relpermdiagnostics.cpp
  tests/test_RestartSerialization.cpp
  tests/test_rftcontainer.cpp
  tests/test_RunningStatistics.cpp
  tests/test_rstconv.cpp
  tests/test_stoppedwells.cpp
  tests/test_timer.cpp
  tests/test_vfpproperties.cpp
  tests/test_wellmodel.cpp
  tests/test_wellprodindexcalculator.cpp
  tests/test_wellstate.cpp
  )

if (HAVE_ECL_INPUT)
  list(APPEND TEST_SOURCE_FILES
    tests/test_nonnc.cpp
    tests/test_GasSatfuncConsistencyChecks.cpp
    tests/test_OilSatfuncConsistencyChecks.cpp
    tests/test_SatfuncCheckPoint.cpp
    tests/test_SatfuncConsistencyCheckManager.cpp
    tests/test_SatfuncConsistencyChecks.cpp
    tests/test_SatfuncConsistencyChecks_parallel.cpp
    tests/test_ThreePointHorizontalSatfuncConsistencyChecks.cpp
    tests/test_WaterSatfuncConsistencyChecks.cpp
  )
endif()

if(MPI_FOUND)
  list(APPEND TEST_SOURCE_FILES tests/test_ghostlastmatrixadapter.cpp
                                tests/test_parallelistlinformation.cpp
                                tests/test_ParallelSerialization.cpp)
endif()

if(CUDA_FOUND)
  if(USE_GPU_BRIDGE)
    list(APPEND TEST_SOURCE_FILES tests/test_cusparseSolver.cpp)
  endif()
endif()

if(CUDA_FOUND OR hip_FOUND)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_converttofloatadapter.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_cublas_handle.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_cublas_safe_call.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuBuffer.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuView.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_cusparse_safe_call.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_gpu_safe_call.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_cuda_check_last_error.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuDILU.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuJac.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuSeqILU0.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_cusparse_handle.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_cuSparse_matrix_operations.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuSparseMatrix.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuSparseTable.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuVector.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_cuVector_operations.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_safe_conversion.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_solver_adapter.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_gpu_ad.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_gpu_linear_two_phase_material.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_gpuPvt.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_gpu_smart_pointers.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_gpu_resources.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_is_gpu_pointer.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_throw_macros_on_gpu.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_blackoilfluidstategpu.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_conditional_storage.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_preconditioner_factory_gpu.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_gpuBlackOilFluidSystem.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuPressureTransferPolicy.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_deviceBlockOperations.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_primaryvarswithdifferentvector.cpp)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_primary_variables_gpu.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_MiniMatrix.cu)
  ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_MiniVector.cu)

  if(MPI_FOUND)
    ADD_CUDA_OR_HIP_FILE(TEST_SOURCE_FILES tests test_GpuOwnerOverlapCopy.cpp)
  endif()

  # for loop providing the flag --expt-relaxed-constexpr to fix some cuda issues with constexpr
  if(NOT CONVERT_CUDA_TO_HIP)
    set(CU_FILES_NEEDING_RELAXED_CONSTEXPR
      tests/gpuistl/test_gpu_ad.cu
      tests/gpuistl/test_gpu_linear_two_phase_material.cu
      tests/gpuistl/test_gpuPvt.cu
      tests/gpuistl/test_gpuBlackOilFluidSystem.cu
      tests/gpuistl/test_GpuSparseMatrix.cu
      tests/gpuistl/test_GpuSparseTable.cu
      tests/gpuistl/test_blackoilfluidstategpu.cu
    )

    foreach(file ${CU_FILES_NEEDING_RELAXED_CONSTEXPR})
        set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "--expt-relaxed-constexpr")
    endforeach()

    set(CU_FILES_NEEDING_FPERMISSIVE
      tests/gpuistl/test_primary_variables_gpu.cu
    )

    foreach(file ${CU_FILES_NEEDING_FPERMISSIVE})
      # Certain structures in OPM requires the -fpermissive flag to compile with nvcc,
      # this enables this for the specific files
      set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "-fpermissive --expt-relaxed-constexpr -Xcompiler=-fpermissive")
    endforeach()
  endif()
endif()

if(USE_GPU_BRIDGE)
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
  tests/norne_pvt_expected.txt
  tests/wells_no_perforation.data
  tests/matr33.txt
  tests/offset_map_matrix.txt
  tests/offset_map_matrix_transposed.txt
  tests/rhs3.txt
  tests/matr33rep.txt
  tests/rhs3rep.txt
  tests/options_flexiblesolver.json
  tests/options_flexiblesolver_simple.json
  tests/GCONSUMP.DATA
  tests/GCONSUMP_COMPLEX.DATA
  tests/GROUP_HIGHER_CONSTRAINTS.DATA
  tests/GLIFT1.DATA
  tests/RC-01_MAST_PRED.DATA
  tests/include/flowl_b_vfp.ecl
  tests/include/flowl_c_vfp.ecl
  tests/include/permx_model5.grdecl
  tests/include/pvt_live_oil_dgas.ecl
  tests/include/relperm.inc
  tests/include/rock.inc
  tests/include/summary.inc
  tests/include/test1_20x30x10.grdecl
  tests/include/well_vfp.ecl
  tests/include/b1_vfp_flowline.inc
  tests/include/d1_vfp_flowline.inc
  tests/include/edit_nnc.inc
  tests/include/flowline_e1_vfp.inc
  tests/include/PVT-WET-GAS.INC
  tests/include/scal_mod2.inc
  tests/include/summary_rc.inc
  tests/test10.partition
  tests/parametersystem.ini
  tests/data/co2injection.dgf
  tests/data/cuvette_11x4.dgf
  tests/data/cuvette_44x24.dgf
  tests/data/fracture.art.dgf
  tests/data/fracture-raw.art
  tests/data/groundwater_1d.dgf
  tests/data/groundwater_2d.dgf
  tests/data/groundwater_3d.dgf
  tests/data/infiltration_50x3.dgf
  tests/data/infiltration_250x20.dgf
  tests/data/obstacle_24x16.dgf
  tests/data/obstacle_48x32.dgf
  tests/data/outflow.dgf
  tests/data/reservoir.dgf
  tests/data/richardslens_24x16.dgf
  tests/data/richardslens_48x32.dgf
  tests/data/richardslens_96x64.dgf
  tests/data/test_stokes.dgf
  tests/data/test_stokes2c.dgf
  tests/data/test_stokes2cni.dgf
  tests/data/waterair.dgf
  )


# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/models/blackoil/blackoilbioeffectsmodules.hh
  opm/models/blackoil/blackoilbioeffectsparams.hpp
  opm/models/blackoil/blackoilboundaryratevector.hh
  opm/models/blackoil/blackoilbrinemodules.hh
  opm/models/blackoil/blackoilbrineparams.hpp
  opm/models/blackoil/blackoilconvectivemixingmodule.hh
  opm/models/blackoil/blackoildarcyfluxmodule.hh
  opm/models/blackoil/blackoildiffusionmodule.hh
  opm/models/blackoil/blackoildispersionmodule.hh
  opm/models/blackoil/blackoilenergymodules.hh
  opm/models/blackoil/blackoilextbomodules.hh
  opm/models/blackoil/blackoilextboparams.hpp
  opm/models/blackoil/blackoilextensivequantities.hh
  opm/models/blackoil/blackoilfoammodules.hh
  opm/models/blackoil/blackoilfoamparams.hpp
  opm/models/blackoil/blackoilvariableandequationindices.hh
  opm/models/blackoil/blackoilintensivequantities.hh
  opm/models/blackoil/blackoillocalresidual.hh
  opm/models/blackoil/blackoillocalresidualtpfa.hh
  opm/models/blackoil/blackoilmeanings.hh
  opm/models/blackoil/blackoilmodel.hh
  opm/models/blackoil/blackoilmoduleparams.hh
  opm/models/blackoil/blackoilnewtonmethod.hpp
  opm/models/blackoil/blackoilnewtonmethodparams.hpp
  opm/models/blackoil/blackoilonephaseindices.hh
  opm/models/blackoil/blackoilpolymermodules.hh
  opm/models/blackoil/blackoilpolymerparams.hpp
  opm/models/blackoil/blackoilprimaryvariables.hh
  opm/models/blackoil/blackoilproblem.hh
  opm/models/blackoil/blackoilproperties.hh
  opm/models/blackoil/blackoilratevector.hh
  opm/models/blackoil/blackoilsolventmodules.hh
  opm/models/blackoil/blackoilsolventparams.hpp
  opm/models/blackoil/blackoiltwophaseindices.hh
  opm/models/common/darcyfluxmodule.hh
  opm/models/common/diffusionmodule.hh
  opm/models/common/directionalmobility.hh
  opm/models/common/energymodule.hh
  opm/models/common/flux.hh
  opm/models/common/forchheimerfluxmodule.hh
  opm/models/common/multiphasebaseextensivequantities.hh
  opm/models/common/multiphasebasemodel.hh
  opm/models/common/multiphasebaseparameters.hh
  opm/models/common/multiphasebaseproblem.hh
  opm/models/common/multiphasebaseproperties.hh
  opm/models/common/quantitycallbacks.hh
  opm/models/common/transfluxmodule.hh
  opm/models/discretefracture/discretefractureextensivequantities.hh
  opm/models/discretefracture/discretefractureintensivequantities.hh
  opm/models/discretefracture/discretefracturelocalresidual.hh
  opm/models/discretefracture/discretefracturemodel.hh
  opm/models/discretefracture/discretefractureprimaryvariables.hh
  opm/models/discretefracture/discretefractureproblem.hh
  opm/models/discretefracture/discretefractureproperties.hh
  opm/models/discretefracture/fracturemapper.hh
  opm/models/discretization/common/baseauxiliarymodule.hh
  opm/models/discretization/common/fvbaseadlocallinearizer.hh
  opm/models/discretization/common/fvbaseboundarycontext.hh
  opm/models/discretization/common/fvbaseconstraints.hh
  opm/models/discretization/common/fvbaseconstraintscontext.hh
  opm/models/discretization/common/fvbasediscretization.hh
  opm/models/discretization/common/fvbasediscretizationfemadapt.hh
  opm/models/discretization/common/fvbaseelementcontext.hh
  opm/models/discretization/common/fvbaseextensivequantities.hh
  opm/models/discretization/common/fvbasefdlocallinearizer.hh
  opm/models/discretization/common/fvbasegradientcalculator.hh
  opm/models/discretization/common/fvbaseintensivequantities.hh
  opm/models/discretization/common/fvbaselinearizer.hh
  opm/models/discretization/common/fvbaselocalresidual.hh
  opm/models/discretization/common/fvbasenewtonconvergencewriter.hh
  opm/models/discretization/common/fvbasenewtonmethod.hh
  opm/models/discretization/common/fvbaseparameters.hh
  opm/models/discretization/common/fvbaseprimaryvariables.hh
  opm/models/discretization/common/fvbaseproblem.hh
  opm/models/discretization/common/fvbaseproperties.hh
  opm/models/discretization/common/linearizationtype.hh
  opm/models/discretization/common/restrictprolong.hh
  opm/models/discretization/common/tpfalinearizer.hh
  opm/models/discretization/ecfv/ecfvbaseoutputmodule.hh
  opm/models/discretization/ecfv/ecfvdiscretization.hh
  opm/models/discretization/ecfv/ecfvgridcommhandlefactory.hh
  opm/models/discretization/ecfv/ecfvproperties.hh
  opm/models/discretization/ecfv/ecfvstencil.hh
  opm/models/discretization/vcfv/vcfvbaseoutputmodule.hh
  opm/models/discretization/vcfv/vcfvdiscretization.hh
  opm/models/discretization/vcfv/vcfvgridcommhandlefactory.hh
  opm/models/discretization/vcfv/vcfvproperties.hh
  opm/models/discretization/vcfv/vcfvstencil.hh
  opm/models/discretization/vcfv/p1fegradientcalculator.hh
  opm/models/flash/flashboundaryratevector.hh
  opm/models/flash/flashextensivequantities.hh
  opm/models/flash/flashindices.hh
  opm/models/flash/flashintensivequantities.hh
  opm/models/flash/flashlocalresidual.hh
  opm/models/flash/flashmodel.hh
  opm/models/flash/flashratevector.hh
  opm/models/flash/flashparameters.hh
  opm/models/flash/flashprimaryvariables.hh
  opm/models/flash/flashproperties.hh
  opm/models/immiscible/immiscibleboundaryratevector.hh
  opm/models/immiscible/immiscibleextensivequantities.hh
  opm/models/immiscible/immiscibleindices.hh
  opm/models/immiscible/immiscibleintensivequantities.hh
  opm/models/immiscible/immisciblelocalresidual.hh
  opm/models/immiscible/immisciblemodel.hh
  opm/models/immiscible/immiscibleprimaryvariables.hh
  opm/models/immiscible/immiscibleproperties.hh
  opm/models/immiscible/immiscibleratevector.hh
  opm/models/io/baseoutputmodule.hh
  opm/models/io/baseoutputwriter.hh
  opm/models/io/basevanguard.hh
  opm/models/io/cubegridvanguard.hh
  opm/models/io/dgfvanguard.hh
  opm/models/io/restart.hpp
  opm/models/io/simplexvanguard.hh
  opm/models/io/structuredgridvanguard.hh
  opm/models/io/unstructuredgridvanguard.hh
  opm/models/io/vtkblackoilbioeffectsmodule.hpp
  opm/models/io/vtkblackoilbioeffectsparams.hpp
  opm/models/io/vtkblackoilenergymodule.hpp
  opm/models/io/vtkblackoilenergyparams.hpp
  opm/models/io/vtkblackoilmodule.hpp
  opm/models/io/vtkblackoilparams.hpp
  opm/models/io/vtkblackoilpolymermodule.hpp
  opm/models/io/vtkblackoilpolymerparams.hpp
  opm/models/io/vtkblackoilsolventmodule.hpp
  opm/models/io/vtkblackoilsolventparams.hpp
  opm/models/io/vtkcompositionmodule.hpp
  opm/models/io/vtkcompositionparams.hpp
  opm/models/io/vtkdiffusionmodule.hpp
  opm/models/io/vtkdiffusionparams.hpp
  opm/models/io/vtkdiscretefracturemodule.hpp
  opm/models/io/vtkdiscretefractureparams.hpp
  opm/models/io/vtkenergymodule.hpp
  opm/models/io/vtkenergyparams.hpp
  opm/models/io/vtkmultiphasemodule.hpp
  opm/models/io/vtkmultiphaseparams.hpp
  opm/models/io/vtkmultiwriter.hh
  opm/models/io/vtkphasepresencemodule.hpp
  opm/models/io/vtkphasepresenceparams.hpp
  opm/models/io/vtkprimaryvarsmodule.hpp
  opm/models/io/vtkprimaryvarsparams.hpp
  opm/models/io/vtkptflashmodule.hpp
  opm/models/io/vtkptflashparams.hpp
  opm/models/io/vtkscalarfunction.hh
  opm/models/io/vtktemperaturemodule.hpp
  opm/models/io/vtktemperatureparams.hpp
  opm/models/io/vtktensorfunction.hh
  opm/models/io/vtkvectorfunction.hh
  opm/models/ncp/ncpboundaryratevector.hh
  opm/models/ncp/ncpextensivequantities.hh
  opm/models/ncp/ncpindices.hh
  opm/models/ncp/ncpintensivequantities.hh
  opm/models/ncp/ncplocalresidual.hh
  opm/models/ncp/ncpmodel.hh
  opm/models/ncp/ncpnewtonmethod.hh
  opm/models/ncp/ncpprimaryvariables.hh
  opm/models/ncp/ncpproperties.hh
  opm/models/ncp/ncpratevector.hh
  opm/models/nonlinear/newtonmethod.hh
  opm/models/nonlinear/newtonmethodparams.hpp
  opm/models/nonlinear/newtonmethodproperties.hh
  opm/models/nonlinear/nullconvergencewriter.hh
  opm/models/parallel/gridcommhandles.hh
  opm/models/parallel/mpibuffer.hh
  opm/models/parallel/tasklets.hpp
  opm/models/parallel/threadedentityiterator.hh
  opm/models/parallel/threadmanager.hpp
  opm/models/ptflash/flashindices.hh
  opm/models/ptflash/flashintensivequantities.hh
  opm/models/ptflash/flashlocalresidual.hh
  opm/models/ptflash/flashmodel.hh
  opm/models/ptflash/flashnewtonmethod.hh
  opm/models/ptflash/flashparameters.hh
  opm/models/ptflash/flashprimaryvariables.hh
  opm/models/pvs/pvsboundaryratevector.hh
  opm/models/pvs/pvsextensivequantities.hh
  opm/models/pvs/pvsindices.hh
  opm/models/pvs/pvsintensivequantities.hh
  opm/models/pvs/pvslocalresidual.hh
  opm/models/pvs/pvsmodel.hh
  opm/models/pvs/pvsnewtonmethod.hh
  opm/models/pvs/pvsprimaryvariables.hh
  opm/models/pvs/pvsproperties.hh
  opm/models/pvs/pvsratevector.hh
  opm/models/richards/richardsboundaryratevector.hh
  opm/models/richards/richardsextensivequantities.hh
  opm/models/richards/richardsindices.hh
  opm/models/richards/richardsintensivequantities.hh
  opm/models/richards/richardslocalresidual.hh
  opm/models/richards/richardsmodel.hh
  opm/models/richards/richardsnewtonmethod.hh
  opm/models/richards/richardsprimaryvariables.hh
  opm/models/richards/richardsproperties.hh
  opm/models/richards/richardsratevector.hh
  opm/models/utils/alignedallocator.hh
  opm/models/utils/basicparameters.hh
  opm/models/utils/basicproperties.hh
  opm/models/utils/genericguard.hh
  opm/models/utils/parametersystem.hpp
  opm/models/utils/pffgridvector.hh
  opm/models/utils/prefetch.hh
  opm/models/utils/propertysystem.hh
  opm/models/utils/quadraturegeometries.hh
  opm/models/utils/signum.hh
  opm/models/utils/simulator.hh
  opm/models/utils/simulatorutils.hpp
  opm/models/utils/start.hh
  opm/models/utils/terminal.hpp
  opm/models/utils/timer.hpp
  opm/models/utils/timerguard.hh
  opm/simulators/flow/ActionHandler.hpp
  opm/simulators/flow/AluGridCartesianIndexMapper.hpp
  opm/simulators/flow/AluGridLevelCartesianIndexMapper.hpp
  opm/simulators/flow/AluGridVanguard.hpp
  opm/simulators/flow/Banners.hpp
  opm/simulators/flow/BaseAquiferModel.hpp
  opm/simulators/flow/BioeffectsContainer.hpp
  opm/simulators/flow/BlackoilModel.hpp
  opm/simulators/flow/BlackoilModel_impl.hpp
  opm/simulators/flow/BlackoilModelConvergenceMonitor.hpp
  opm/simulators/flow/BlackoilModelNldd.hpp
  opm/simulators/flow/BlackoilModelParameters.hpp
  opm/simulators/flow/BlackoilModelProperties.hpp
  opm/simulators/flow/CO2H2Container.hpp
  opm/simulators/flow/CollectDataOnIORank.hpp
  opm/simulators/flow/CollectDataOnIORank_impl.hpp
  opm/simulators/flow/CompositionalContainer.hpp
  opm/simulators/flow/ConvergenceOutputConfiguration.hpp
  opm/simulators/flow/countGlobalCells.hpp
  opm/simulators/flow/CpGridVanguard.hpp
  opm/simulators/flow/DummyGradientCalculator.hpp
  opm/simulators/flow/EclGenericWriter.hpp
  opm/simulators/flow/EclGenericWriter_impl.hpp
  opm/simulators/flow/EclWriter.hpp
  opm/simulators/flow/EquilInitializer.hpp
  opm/simulators/flow/ExtboContainer.hpp
  opm/simulators/flow/ExtraConvergenceOutputThread.hpp
  opm/simulators/flow/FemCpGridCompat.hpp
  opm/simulators/flow/FIBlackoilModel.hpp
  opm/simulators/flow/FIPContainer.hpp
  opm/simulators/flow/FlowBaseProblemProperties.hpp
  opm/simulators/flow/FlowBaseVanguard.hpp
  opm/simulators/flow/FlowGenericProblem.hpp
  opm/simulators/flow/FlowGenericProblem_impl.hpp
  opm/simulators/flow/FlowGenericVanguard.hpp
  opm/simulators/flow/FlowMain.hpp
  opm/simulators/flow/FlowProblem.hpp
  opm/simulators/flow/FlowProblemBlackoil.hpp
  opm/simulators/flow/FlowProblemBlackoilProperties.hpp
  opm/simulators/flow/FlowProblemComp.hpp
  opm/simulators/flow/FlowProblemCompProperties.hpp
  opm/simulators/flow/FlowProblemParameters.hpp
  opm/simulators/flow/FlowsContainer.hpp
  opm/simulators/flow/FlowUtils.hpp
  opm/simulators/flow/FlowsData.hpp
  opm/simulators/flow/FlowThresholdPressure.hpp
  opm/simulators/flow/GenericCpGridVanguard.hpp
  opm/simulators/flow/GenericOutputBlackoilModule.hpp
  opm/simulators/flow/GenericTemperatureModel.hpp
  opm/simulators/flow/GenericTemperatureModel_impl.hpp
  opm/simulators/flow/GenericThresholdPressure.hpp
  opm/simulators/flow/GenericThresholdPressure_impl.hpp
  opm/simulators/flow/GenericTracerModel.hpp
  opm/simulators/flow/GenericTracerModel_impl.hpp
  opm/simulators/flow/HybridNewton.hpp
  opm/simulators/flow/HybridNewtonConfig.hpp
  opm/simulators/flow/InterRegFlows.hpp
  opm/simulators/flow/KeywordValidation.hpp
  opm/simulators/flow/LogOutputHelper.hpp
  opm/simulators/flow/Main.hpp
  opm/simulators/flow/MechContainer.hpp
  opm/simulators/flow/MixingRateControls.hpp
  opm/simulators/flow/NewTranFluxModule.hpp
  opm/simulators/flow/NlddReporting.hpp
  opm/simulators/flow/NonlinearSolver.hpp
  opm/simulators/flow/OutputBlackoilModule.hpp
  opm/simulators/flow/OutputCompositionalModule.hpp
  opm/simulators/flow/OutputExtractor.hpp
  opm/simulators/flow/partitionCells.hpp
  opm/simulators/flow/PolyhedralGridVanguard.hpp
  opm/simulators/flow/priVarsPacking.hpp
  opm/simulators/flow/RFTContainer.hpp
  opm/simulators/flow/RSTConv.hpp
  opm/simulators/flow/RegionPhasePVAverage.hpp
  opm/simulators/flow/SimulatorConvergenceOutput.hpp
  opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp
  opm/simulators/flow/SimulatorReportBanners.hpp
  opm/simulators/flow/SimulatorSerializer.hpp
  opm/simulators/flow/SolutionContainers.hpp
  opm/simulators/flow/SubDomain.hpp
  opm/simulators/flow/TTagFlowProblemTPFA.hpp
  opm/simulators/flow/TTagFlowProblemGasWater.hpp
  opm/simulators/flow/TTagFlowProblemOnePhase.hpp
  opm/simulators/flow/TracerContainer.hpp
  opm/simulators/flow/TemperatureModel.hpp
  opm/simulators/flow/TracerModel.hpp
  opm/simulators/flow/Transmissibility.hpp
  opm/simulators/flow/Transmissibility_impl.hpp
  opm/simulators/flow/ValidationFunctions.hpp
  opm/simulators/flow/VtkTracerModule.hpp
  opm/simulators/flow/equil/EquilibrationHelpers.hpp
  opm/simulators/flow/equil/EquilibrationHelpers_impl.hpp
  opm/simulators/flow/equil/InitStateEquil.hpp
  opm/simulators/flow/equil/InitStateEquil_impl.hpp
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
  opm/simulators/aquifers/SupportsFaceTag.hpp
  opm/simulators/linalg/AbstractISTLSolver.hpp
  opm/simulators/linalg/amgcpr.hh
  opm/simulators/linalg/bicgstabsolver.hh
  opm/simulators/linalg/blacklist.hh
  opm/simulators/linalg/combinedcriterion.hh
  opm/simulators/linalg/convergencecriterion.hh
  opm/simulators/linalg/DILU.hpp
  opm/simulators/linalg/domesticoverlapfrombcrsmatrix.hh
  opm/simulators/linalg/elementborderlistfromgrid.hh
  opm/simulators/linalg/exportSystem.hpp
  opm/simulators/linalg/extractMatrix.hpp
  opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp
  opm/simulators/linalg/ExtraSmoothers.hpp
  opm/simulators/linalg/findOverlapRowsAndColumns.hpp
  opm/simulators/linalg/fixpointcriterion.hh
  opm/simulators/linalg/FlexibleSolver.hpp
  opm/simulators/linalg/FlexibleSolver_impl.hpp
  opm/simulators/linalg/FlowLinearSolverParameters.hpp
  opm/simulators/linalg/foreignoverlapfrombcrsmatrix.hh
  opm/simulators/linalg/getQuasiImpesWeights.hpp
  opm/simulators/linalg/globalindices.hh
  opm/simulators/linalg/GraphColoring.hpp
  opm/simulators/linalg/ilufirstelement.hh
  opm/simulators/linalg/is_gpu_operator.hpp
  opm/simulators/linalg/ISTLSolver.hpp
  opm/simulators/linalg/ISTLSolverRuntimeOptionProxy.hpp
  opm/simulators/linalg/istlpreconditionerwrappers.hh
  opm/simulators/linalg/istlsolverwrappers.hh
  opm/simulators/linalg/istlsparsematrixadapter.hh
  opm/simulators/linalg/linalgparameters.hh
  opm/simulators/linalg/linalgproperties.hh
  opm/simulators/linalg/LinearSolverAcceleratorType.hpp
  opm/simulators/linalg/linearsolverreport.hh
  opm/simulators/linalg/matrixblock.hh
  opm/simulators/linalg/MatrixMarketSpecializations.hpp
  opm/simulators/linalg/nullborderlistmanager.hh
  opm/simulators/linalg/overlappingbcrsmatrix.hh
  opm/simulators/linalg/overlappingblockvector.hh
  opm/simulators/linalg/overlappingoperator.hh
  opm/simulators/linalg/overlappingpreconditioner.hh
  opm/simulators/linalg/overlappingscalarproduct.hh
  opm/simulators/linalg/overlaptypes.hh
  opm/simulators/linalg/OwningBlockPreconditioner.hpp
  opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp
  opm/simulators/linalg/MILU.hpp
  opm/simulators/linalg/parallelamgbackend.hh
  opm/simulators/linalg/parallelbasebackend.hh
  opm/simulators/linalg/parallelbicgstabbackend.hh
  opm/simulators/linalg/parallelistlbackend.hh
  opm/simulators/linalg/ParallelIstlInformation.hpp
  opm/simulators/linalg/ParallelOverlappingILU0.hpp
  opm/simulators/linalg/ParallelRestrictedAdditiveSchwarz.hpp
  opm/simulators/linalg/PreconditionerFactoryGPUIncludeWrapper.hpp
  opm/simulators/linalg/PreconditionerFactory.hpp
  opm/simulators/linalg/PreconditionerFactory_impl.hpp
  opm/simulators/linalg/printlinearsolverparameter.hpp
  opm/simulators/linalg/StandardPreconditioners.hpp
  opm/simulators/linalg/StandardPreconditioners_mpi.hpp
  opm/simulators/linalg/StandardPreconditioners_serial.hpp
  opm/simulators/linalg/StandardPreconditioners_gpu_serial.hpp
  opm/simulators/linalg/StandardPreconditioners_gpu_mpi.hpp
  opm/simulators/linalg/PreconditionerWithUpdate.hpp
  opm/simulators/linalg/PressureBhpTransferPolicy.hpp
  opm/simulators/linalg/PressureSolverPolicy.hpp
  opm/simulators/linalg/PressureTransferPolicy.hpp
  opm/simulators/linalg/PropertyTree.hpp
  opm/simulators/linalg/residreductioncriterion.hh
  opm/simulators/linalg/SmallDenseMatrixUtils.hpp
  opm/simulators/linalg/setupPropertyTree.hpp
  opm/simulators/linalg/superlubackend.hh
  opm/simulators/linalg/twolevelmethodcpr.hh
  opm/simulators/linalg/vertexborderlistfromgrid.hh
  opm/simulators/linalg/weightedresidreductioncriterion.hh
  opm/simulators/linalg/WellOperators.hpp
  opm/simulators/linalg/WriteSystemMatrixHelper.hpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp
  opm/simulators/timestepping/AdaptiveTimeStepping.hpp
  opm/simulators/timestepping/AdaptiveTimeStepping_impl.hpp
  opm/simulators/timestepping/ConvergenceReport.hpp
  opm/simulators/timestepping/EclTimeSteppingParams.hpp
  opm/simulators/timestepping/TimeStepControl.hpp
  opm/simulators/timestepping/TimeStepControlInterface.hpp
  opm/simulators/timestepping/SimulatorTimer.hpp
  opm/simulators/timestepping/SimulatorReport.hpp
  opm/simulators/timestepping/SimulatorTimerInterface.hpp
  opm/simulators/timestepping/gatherConvergenceReport.hpp
  opm/simulators/utils/ComponentName.hpp
  opm/simulators/utils/DeferredLogger.hpp
  opm/simulators/utils/DeferredLoggingErrorHelpers.hpp
  opm/simulators/utils/ParallelEclipseState.hpp
  opm/simulators/utils/ParallelFileMerger.hpp
  opm/simulators/utils/ParallelNLDDPartitioningZoltan.hpp
  opm/simulators/utils/ParallelRestart.hpp
  opm/simulators/utils/PressureAverage.hpp
  opm/simulators/utils/PropsDataHandle.hpp
  opm/simulators/utils/SerializationPackers.hpp
  opm/simulators/utils/VectorVectorDataHandle.hpp
  opm/simulators/utils/compressPartition.hpp
  opm/simulators/utils/gatherDeferredLogger.hpp
  opm/simulators/utils/moduleVersion.hpp
  opm/simulators/utils/ParallelCommunication.hpp
  opm/simulators/utils/ParallelSerialization.hpp
  opm/simulators/utils/readDeck.hpp
  opm/simulators/utils/satfunc/RelpermDiagnostics.hpp
  opm/simulators/wells/ALQState.hpp
  opm/simulators/wells/BlackoilWellModel.hpp
  opm/simulators/wells/BlackoilWellModel_impl.hpp
  opm/simulators/wells/BlackoilWellModelConstraints.hpp
  opm/simulators/wells/BlackoilWellModelGasLift.hpp
  opm/simulators/wells/BlackoilWellModelGasLift_impl.hpp
  opm/simulators/wells/BlackoilWellModelGeneric.hpp
  opm/simulators/wells/BlackoilWellModelGuideRates.hpp
  opm/simulators/wells/BlackoilWellModelNetwork.hpp
  opm/simulators/wells/BlackoilWellModelNetwork_impl.hpp
  opm/simulators/wells/BlackoilWellModelNetworkGeneric.hpp
  opm/simulators/wells/BlackoilWellModelNldd.hpp
  opm/simulators/wells/BlackoilWellModelNldd_impl.hpp
  opm/simulators/wells/BlackoilWellModelRestart.hpp
  opm/simulators/wells/BlackoilWellModelWBP.hpp
  opm/simulators/wells/ConnectionIndexMap.hpp
  opm/simulators/wells/ConnFiltrateData.hpp
  opm/simulators/wells/ConnFracStatistics.hpp
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
  opm/simulators/wells/GroupStateHelper.hpp
  opm/simulators/wells/GuideRateHandler.hpp
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
  opm/simulators/wells/RatioCalculator.hpp
  opm/simulators/wells/RegionAttributeHelpers.hpp
  opm/simulators/wells/RegionAverageCalculator.hpp
  opm/simulators/wells/RunningStatistics.hpp
  opm/simulators/wells/RuntimePerforation.hpp
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
  opm/simulators/wells/WellHelpers.hpp
  opm/simulators/wells/WellInterfaceFluidSystem.hpp
  opm/simulators/wells/WellInterfaceGeneric.hpp
  opm/simulators/wells/WellInterface.hpp
  opm/simulators/wells/WellInterface_impl.hpp
  opm/simulators/wells/WellInterfaceIndices.hpp
  opm/simulators/wells/WellProdIndexCalculator.hpp
  opm/simulators/wells/WellState.hpp
  opm/simulators/wells/WellTest.hpp
  opm/simulators/wells/WellTracerRate.hpp
  opm/simulators/wells/WGState.hpp
)
if (USE_GPU_BRIDGE)
  list (APPEND PUBLIC_HEADER_FILES
    opm/simulators/linalg/gpubridge/amgclSolverBackend.hpp
    opm/simulators/linalg/gpubridge/GpuBridge.hpp
    opm/simulators/linalg/gpubridge/GpuResult.hpp
    opm/simulators/linalg/gpubridge/GpuSolver.hpp
    opm/simulators/linalg/gpubridge/CprCreation.hpp
    opm/simulators/linalg/gpubridge/Preconditioner.hpp
    opm/simulators/linalg/gpubridge/Misc.hpp
    opm/simulators/linalg/gpubridge/opencl/openclBILU0.hpp
    opm/simulators/linalg/gpubridge/BlockedMatrix.hpp
    opm/simulators/linalg/gpubridge/opencl/openclCPR.hpp
    opm/simulators/linalg/gpubridge/cuda/cuda_header.hpp
    opm/simulators/linalg/gpubridge/cuda/cusparseSolverBackend.hpp
    opm/simulators/linalg/gpubridge/opencl/ChowPatelIlu.hpp
    opm/simulators/linalg/gpubridge/opencl/openclBISAI.hpp
    opm/simulators/linalg/gpubridge/Reorder.hpp
    opm/simulators/linalg/gpubridge/opencl/opencl.hpp
    opm/simulators/linalg/gpubridge/opencl/openclKernels.hpp
    opm/simulators/linalg/gpubridge/opencl/OpenclMatrix.hpp
    opm/simulators/linalg/gpubridge/opencl/openclPreconditioner.hpp
    opm/simulators/linalg/gpubridge/opencl/openclSolverBackend.hpp
    opm/simulators/linalg/gpubridge/opencl/openclWellContributions.hpp
    opm/simulators/linalg/gpubridge/Matrix.hpp
    opm/simulators/linalg/gpubridge/MultisegmentWellContribution.hpp
    opm/simulators/linalg/gpubridge/rocm/hipKernels.hpp
    opm/simulators/linalg/gpubridge/rocm/rocalutionSolverBackend.hpp
    opm/simulators/linalg/gpubridge/rocm/rocsparseBILU0.hpp
    opm/simulators/linalg/gpubridge/rocm/rocsparseCPR.hpp
    opm/simulators/linalg/gpubridge/rocm/rocsparsePreconditioner.hpp
    opm/simulators/linalg/gpubridge/rocm/rocsparseSolverBackend.hpp
    opm/simulators/linalg/gpubridge/rocm/rocsparseWellContributions.hpp
    opm/simulators/linalg/gpubridge/rocm/rocsparseMatrix.hpp
    opm/simulators/linalg/gpubridge/WellContributions.hpp
    opm/simulators/linalg/ISTLSolverGpuBridge.hpp
  )
endif()

if (HAVE_ECL_INPUT)
  list (APPEND PUBLIC_HEADER_FILES
    opm/simulators/utils/satfunc/GasPhaseConsistencyChecks.hpp
    opm/simulators/utils/satfunc/OilPhaseConsistencyChecks.hpp
    opm/simulators/utils/satfunc/PhaseCheckBase.hpp
    opm/simulators/utils/satfunc/SatfuncCheckPointInterface.hpp
    opm/simulators/utils/satfunc/SatfuncConsistencyCheckManager.hpp
    opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp
    opm/simulators/utils/satfunc/ScaledSatfuncCheckPoint.hpp
    opm/simulators/utils/satfunc/ThreePointHorizontalConsistencyChecks.hpp
    opm/simulators/utils/satfunc/UnscaledSatfuncCheckPoint.hpp
    opm/simulators/utils/satfunc/WaterPhaseConsistencyChecks.hpp
  )
endif()

if (HAVE_AVX2_EXTENSION)
  list (APPEND PUBLIC_HEADER_FILES
    opm/simulators/linalg/mixed/bslv.h
    opm/simulators/linalg/mixed/bsr.h
    opm/simulators/linalg/mixed/prec.h
    opm/simulators/linalg/mixed/vec.h
    opm/simulators/linalg/mixed/wrapper.hpp
  )
endif()

if (Damaris_FOUND AND MPI_FOUND AND USE_DAMARIS_LIB)
  list (APPEND PUBLIC_HEADER_FILES
    opm/simulators/utils/DamarisKeywords.hpp
    opm/simulators/utils/DamarisOutputModule.hpp
    opm/simulators/flow/DamarisParameters.hpp
    opm/simulators/flow/DamarisWriter.hpp
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
  examples/art2dgf.cpp
  examples/co2injection_flash_ecfv.cpp
  examples/co2injection_flash_ni_ecfv.cpp
  examples/co2injection_flash_ni_vcfv.cpp
  examples/co2injection_flash_vcfv.cpp
  examples/co2injection_immiscible_ecfv.cpp
  examples/co2injection_immiscible_ni_ecfv.cpp
  examples/co2injection_immiscible_ni_vcfv.cpp
  examples/co2injection_immiscible_vcfv.cpp
  examples/co2injection_ncp_ecfv.cpp
  examples/co2injection_ncp_ni_ecfv.cpp
  examples/co2injection_ncp_ni_vcfv.cpp
  examples/co2injection_ncp_vcfv.cpp
  examples/co2injection_pvs_ecfv.cpp
  examples/co2injection_pvs_ni_ecfv.cpp
  examples/co2injection_pvs_ni_vcfv.cpp
  examples/co2_ptflash_ecfv.cpp
  examples/co2injection_pvs_vcfv.cpp
  examples/cuvette_pvs.cpp
  examples/diffusion_flash.cpp
  examples/diffusion_ncp.cpp
  examples/diffusion_pvs.cpp
  examples/groundwater_immiscible.cpp
  examples/infiltration_pvs.cpp
  examples/lens_immiscible_ecfv_ad.cpp
  examples/lens_immiscible_ecfv_ad_23.cpp
  examples/lens_immiscible_ecfv_ad_trans.cpp
  examples/lens_immiscible_vcfv_ad.cpp
  examples/lens_immiscible_vcfv_fd.cpp
  examples/lens_richards_ecfv.cpp
  examples/lens_richards_vcfv.cpp
  examples/obstacle_immiscible.cpp
  examples/obstacle_ncp.cpp
  examples/obstacle_pvs.cpp
  examples/outflow_pvs.cpp
  examples/powerinjection_darcy_ad.cpp
  examples/powerinjection_darcy_fd.cpp
  examples/powerinjection_forchheimer_ad.cpp
  examples/powerinjection_forchheimer_fd.cpp
  examples/reservoir_blackoil_ecfv.cpp
  examples/reservoir_blackoil_vcfv.cpp
  examples/reservoir_ncp_ecfv.cpp
  examples/reservoir_ncp_vcfv.cpp
  examples/printvfp.cpp
  examples/tutorial1.cpp
  examples/waterair_pvs_ni.cpp
)
if(HDF5_FOUND)
  list (APPEND EXAMPLE_SOURCE_FILES
    examples/opmrst_inspect.cpp
  )
endif()
if(dune-alugrid_FOUND)
  list (APPEND EXAMPLE_SOURCE_FILES
    examples/finger_immiscible_ecfv.cpp
    examples/finger_immiscible_vcfv.cpp
    examples/fracture_discretefracture.cpp
  )
endif()
if(MPI_FOUND)
  list (APPEND MAIN_SOURCE_FILES
    opm/simulators/flow/rescoup/ReservoirCoupling.cpp
    opm/simulators/flow/rescoup/ReservoirCouplingMaster.cpp
    opm/simulators/flow/rescoup/ReservoirCouplingMasterReportStep.cpp
    opm/simulators/flow/rescoup/ReservoirCouplingSlave.cpp
    opm/simulators/flow/rescoup/ReservoirCouplingSlaveReportStep.cpp
    opm/simulators/flow/rescoup/ReservoirCouplingSpawnSlaves.cpp
    opm/simulators/flow/rescoup/ReservoirCouplingTimeStepper.cpp
    opm/simulators/wells/GroupTargetCalculator.cpp
    opm/simulators/wells/rescoup/RescoupReceiveGroupTargets.cpp
    opm/simulators/wells/rescoup/RescoupReceiveSlaveGroupData.cpp
    opm/simulators/wells/rescoup/RescoupSendSlaveGroupData.cpp
    opm/simulators/wells/rescoup/RescoupTargetCalculator.cpp
  )
  list (APPEND PUBLIC_HEADER_FILES
    opm/simulators/flow/rescoup/ReservoirCoupling.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingEnabled.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingErrorMacros.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingMasterReportStep.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingSlaveReportStep.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingSpawnSlaves.hpp
    opm/simulators/flow/rescoup/ReservoirCouplingTimeStepper.hpp
    opm/simulators/wells/GroupTargetCalculator.hpp
    opm/simulators/wells/rescoup/RescoupProxy.hpp
    opm/simulators/wells/rescoup/RescoupReceiveSlaveGroupData.hpp
    opm/simulators/wells/rescoup/RescoupReceiveGroupTargets.hpp
    opm/simulators/wells/rescoup/RescoupSendSlaveGroupData.hpp
    opm/simulators/wells/rescoup/RescoupTargetCalculator.hpp
    )
  list (APPEND TEST_SOURCE_FILES
    tests/rescoup/test_chopstep.cpp
  )
endif()
if(HYPRE_FOUND)
  list(APPEND PUBLIC_HEADER_FILES
    opm/simulators/linalg/HyprePreconditioner.hpp
  )
endif()

if(AMGX_FOUND)
  list(APPEND PUBLIC_HEADER_FILES
    opm/simulators/linalg/AmgxPreconditioner.hpp
  )
endif()

if (CONVERT_CUDA_TO_HIP)
  add_custom_target(hipified_headers  DEPENDS ${PUBLIC_HEADER_FILES_HIPIFIED})
endif()
