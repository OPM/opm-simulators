# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

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
  opm/autodiff/NewtonIterationBlackoilCPR.cpp
  opm/autodiff/NewtonIterationBlackoilInterleaved.cpp
  opm/autodiff/NewtonIterationBlackoilSimple.cpp
  opm/autodiff/NewtonIterationUtilities.cpp
  opm/autodiff/GridHelpers.cpp
  opm/autodiff/ImpesTPFAAD.cpp
  opm/autodiff/moduleVersion.cpp
  opm/autodiff/multiPhaseUpwind.cpp
  opm/autodiff/SimulatorFullyImplicitBlackoilOutput.cpp
  opm/autodiff/SimulatorIncompTwophaseAd.cpp
  opm/autodiff/TransportSolverTwophaseAd.cpp
  opm/autodiff/BlackoilPropsAdFromDeck.cpp
  opm/autodiff/SolventPropsAdFromDeck.cpp
  opm/autodiff/BlackoilModelParameters.cpp
  opm/autodiff/WellDensitySegmented.cpp
  opm/autodiff/LinearisedBlackoilResidual.cpp
  opm/autodiff/VFPProperties.cpp
  opm/autodiff/VFPProdProperties.cpp
  opm/autodiff/VFPInjProperties.cpp
  opm/autodiff/WellMultiSegment.cpp
  opm/autodiff/MultisegmentWells.cpp
  opm/autodiff/MissingFeatures.cpp
  opm/polymer/PolymerState.cpp
  opm/polymer/PolymerBlackoilState.cpp
  opm/polymer/CompressibleTpfaPolymer.cpp
  opm/polymer/IncompTpfaPolymer.cpp
  opm/polymer/PolymerInflow.cpp
  opm/polymer/PolymerProperties.cpp
  opm/polymer/polymerUtilities.cpp
  opm/polymer/SimulatorCompressiblePolymer.cpp
  opm/polymer/SimulatorPolymer.cpp
  opm/polymer/TransportSolverTwophaseCompressiblePolymer.cpp
  opm/polymer/TransportSolverTwophasePolymer.cpp
  opm/polymer/fullyimplicit/PolymerPropsAd.cpp
	opm/simulators/flow_ebos_blackoil.cpp
  opm/simulators/flow_ebos_gasoil.cpp
  opm/simulators/flow_ebos_oilwater.cpp
  opm/simulators/flow_ebos_polymer.cpp
  opm/simulators/flow_ebos_solvent.cpp
  opm/simulators/ensureDirectoryExists.cpp
  opm/simulators/SimulatorCompressibleTwophase.cpp
  opm/simulators/SimulatorIncompTwophase.cpp
  opm/simulators/WellSwitchingLogger.cpp
  opm/simulators/vtk/writeVtkData.cpp
  opm/simulators/timestepping/TimeStepControl.cpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.cpp
  opm/simulators/timestepping/SimulatorTimer.cpp
  )


# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/test_autodiffhelpers.cpp
  tests/test_autodiffmatrix.cpp
  tests/test_block.cpp
  tests/test_boprops_ad.cpp
  tests/test_rateconverter.cpp
  tests/test_span.cpp
  tests/test_syntax.cpp
  tests/test_scalar_mult.cpp
  tests/test_transmissibilitymultipliers.cpp
  tests/test_welldensitysegmented.cpp
  tests/test_vfpproperties.cpp
  tests/test_singlecellsolves.cpp
  tests/test_solventprops_ad.cpp
  tests/test_multisegmentwells.cpp
  tests/test_multiphaseupwind.cpp
  tests/test_wellmodel.cpp
  # tests/test_thresholdpressure.cpp
  tests/test_wellswitchlogger.cpp
  tests/test_timer.cpp
  tests/test_invert.cpp
  )

list (APPEND TEST_DATA_FILES
  tests/fluid.data
  tests/VFPPROD1
  tests/VFPPROD2
  tests/msw.data
  tests/TESTTIMER.DATA
  tests/TESTWELLMODEL.DATA
  )


# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
  examples/find_zero.cpp
  examples/flow_legacy.cpp
  examples/flow_reorder.cpp
  examples/flow_sequential.cpp
  examples/flow.cpp
  examples/flow_ebos.cpp
  examples/flow_ebos_2p.cpp
  examples/flow_ebos_solvent.cpp
  examples/flow_ebos_polymer.cpp
  examples/flow_multisegment.cpp
  examples/flow_solvent.cpp
  examples/sim_2p_incomp.cpp
  examples/sim_2p_incomp_ad.cpp
  examples/sim_2p_comp_reorder.cpp
  examples/sim_simple.cpp
  examples/opm_init_check.cpp
  examples/sim_poly2p_comp_reorder.cpp
  examples/sim_poly2p_incomp_reorder.cpp
  examples/flow_polymer.cpp
  examples/wells_example.cpp
  )

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
  examples/sim_2p_incomp.cpp
  examples/sim_2p_incomp_ad.cpp
  examples/sim_2p_comp_reorder.cpp
  examples/flow.cpp
  examples/flow_ebos.cpp
  examples/flow_ebos_2p.cpp
  examples/flow_ebos_solvent.cpp
  examples/flow_ebos_polymer.cpp
  examples/flow_legacy.cpp
  examples/flow_reorder.cpp
  examples/flow_sequential.cpp
  examples/flow_solvent.cpp
  examples/opm_init_check.cpp
  examples/sim_poly2p_comp_reorder.cpp
  examples/sim_poly2p_incomp_reorder.cpp
  examples/flow_polymer.cpp
  )

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/autodiff/AdditionalObjectDeleter.hpp
  opm/autodiff/AutoDiffBlock.hpp
  opm/autodiff/AutoDiffHelpers.hpp
  opm/autodiff/AutoDiffMatrix.hpp
  opm/autodiff/AutoDiff.hpp
  opm/autodiff/BackupRestore.hpp
  opm/autodiff/BlackoilDetails.hpp
  opm/autodiff/BlackoilLegacyDetails.hpp
  opm/autodiff/BlackoilModel.hpp
  opm/autodiff/BlackoilModelBase.hpp
  opm/autodiff/BlackoilModelBase_impl.hpp
  opm/autodiff/BlackoilModelEnums.hpp
  opm/autodiff/BlackoilModelParameters.hpp
  opm/autodiff/BlackoilPressureModel.hpp
  opm/autodiff/BlackoilPropsAdFromDeck.hpp
  opm/autodiff/SolventPropsAdFromDeck.hpp
  opm/autodiff/Compat.hpp
  opm/autodiff/CPRPreconditioner.hpp
  opm/autodiff/createGlobalCellArray.hpp
  opm/autodiff/DefaultBlackoilSolutionState.hpp
  opm/autodiff/BlackoilSequentialModel.hpp
  opm/autodiff/BlackoilSolventModel.hpp
  opm/autodiff/BlackoilSolventModel_impl.hpp
  opm/autodiff/BlackoilMultiSegmentModel.hpp
  opm/autodiff/BlackoilMultiSegmentModel_impl.hpp
  opm/autodiff/BlackoilReorderingTransportModel.hpp
  opm/autodiff/BlackoilTransportModel.hpp
  opm/autodiff/fastSparseOperations.hpp
  opm/autodiff/DebugTimeReport.hpp
  opm/autodiff/DuneMatrix.hpp
  opm/autodiff/ExtractParallelGridInformationToISTL.hpp
  opm/autodiff/FlowMain.hpp
  opm/autodiff/FlowMainEbos.hpp
  opm/autodiff/FlowMainPolymer.hpp
  opm/autodiff/FlowMainSequential.hpp
  opm/autodiff/FlowMainSolvent.hpp
  opm/autodiff/GeoProps.hpp
  opm/autodiff/GridHelpers.hpp
  opm/autodiff/GridInit.hpp
  opm/autodiff/ImpesTPFAAD.hpp
  opm/autodiff/ISTLSolver.hpp
  opm/autodiff/IterationReport.hpp
  opm/autodiff/moduleVersion.hpp
  opm/autodiff/multiPhaseUpwind.hpp
  opm/autodiff/NewtonIterationBlackoilCPR.hpp
  opm/autodiff/NewtonIterationBlackoilInterface.hpp
  opm/autodiff/NewtonIterationBlackoilInterleaved.hpp
  opm/autodiff/NewtonIterationBlackoilSimple.hpp
  opm/autodiff/NewtonIterationUtilities.hpp
  opm/autodiff/NonlinearSolver.hpp
  opm/autodiff/NonlinearSolver_impl.hpp
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
  opm/autodiff/SimulatorFullyImplicitBlackoilSolvent.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilSolvent_impl.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilMultiSegment.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilMultiSegment_impl.hpp
  opm/autodiff/SimulatorIncompTwophaseAd.hpp
  opm/autodiff/SimulatorSequentialBlackoil.hpp
  opm/autodiff/TransportSolverTwophaseAd.hpp
  opm/autodiff/WellDensitySegmented.hpp
  opm/autodiff/WellStateFullyImplicitBlackoil.hpp
  opm/autodiff/WellStateFullyImplicitBlackoilSolvent.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp
  opm/autodiff/VFPProperties.hpp
  opm/autodiff/VFPHelpers.hpp
  opm/autodiff/VFPProdProperties.hpp
  opm/autodiff/VFPInjProperties.hpp
  opm/autodiff/WellStateMultiSegment.hpp
  opm/autodiff/WellMultiSegment.hpp
  opm/autodiff/MultisegmentWells.hpp
  opm/autodiff/MultisegmentWells_impl.hpp
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
  opm/autodiff/StandardWellsSolvent.hpp
  opm/autodiff/StandardWellsSolvent_impl.hpp
  opm/autodiff/MissingFeatures.hpp
  opm/autodiff/ThreadHandle.hpp
  opm/polymer/CompressibleTpfaPolymer.hpp
  opm/polymer/GravityColumnSolverPolymer.hpp
  opm/polymer/GravityColumnSolverPolymer_impl.hpp
  opm/polymer/IncompPropertiesDefaultPolymer.hpp
  opm/polymer/IncompTpfaPolymer.hpp
  opm/polymer/PolymerBlackoilState.hpp
  opm/polymer/PolymerInflow.hpp
  opm/polymer/PolymerProperties.hpp
  opm/polymer/PolymerState.hpp
  opm/polymer/polymerUtilities.hpp
  opm/polymer/SimulatorCompressiblePolymer.hpp
  opm/polymer/SimulatorPolymer.hpp
  opm/polymer/SinglePointUpwindTwoPhasePolymer.hpp
  opm/polymer/TransportSolverTwophaseCompressiblePolymer.hpp
  opm/polymer/Point2D.hpp
  opm/polymer/TransportSolverTwophasePolymer.hpp
  opm/polymer/fullyimplicit/PolymerPropsAd.hpp
  opm/polymer/fullyimplicit/BlackoilPolymerModel.hpp
  opm/polymer/fullyimplicit/BlackoilPolymerModel_impl.hpp
  opm/polymer/fullyimplicit/SimulatorFullyImplicitBlackoilPolymer.hpp
  opm/polymer/fullyimplicit/SimulatorFullyImplicitBlackoilPolymer_impl.hpp
  opm/polymer/fullyimplicit/WellStateFullyImplicitBlackoilPolymer.hpp
	opm/simulators/flow_ebos_blackoil.hpp
  opm/simulators/flow_ebos_gasoil.hpp
  opm/simulators/flow_ebos_oilwater.hpp
  opm/simulators/flow_ebos_polymer.hpp
  opm/simulators/flow_ebos_solvent.hpp
  opm/simulators/ensureDirectoryExists.hpp
  opm/simulators/ParallelFileMerger.hpp
  opm/simulators/SimulatorCompressibleTwophase.hpp
  opm/simulators/SimulatorIncompTwophase.hpp
  opm/simulators/thresholdPressures.hpp
  opm/simulators/WellSwitchingLogger.hpp
  opm/simulators/vtk/writeVtkData.hpp
  opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp
  opm/simulators/timestepping/AdaptiveTimeStepping.hpp
  opm/simulators/timestepping/AdaptiveTimeStepping_impl.hpp
  opm/simulators/timestepping/TimeStepControl.hpp
  opm/simulators/timestepping/TimeStepControlInterface.hpp
  opm/simulators/timestepping/SimulatorTimer.hpp
  opm/simulators/timestepping/SimulatorTimerInterface.hpp
  )

