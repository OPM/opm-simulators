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
	# place the flow_ebos_*.cpp files on top of the list because they
	# take the longest to compile, and compiling them first speeds up
	# parallel builds because it allows the jobserver to do better scheduling
  opm/simulators/flow_ebos_blackoil.cpp
  opm/simulators/flow_ebos_gasoil.cpp
  opm/simulators/flow_ebos_oilwater.cpp
  opm/simulators/flow_ebos_polymer.cpp
  opm/simulators/flow_ebos_solvent.cpp

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
  opm/autodiff/BlackoilModelParameters.cpp
  opm/autodiff/WellDensitySegmented.cpp
  opm/autodiff/LinearisedBlackoilResidual.cpp
  opm/autodiff/VFPProperties.cpp
  opm/autodiff/VFPProdProperties.cpp
  opm/autodiff/VFPInjProperties.cpp
  opm/autodiff/MissingFeatures.cpp
  opm/core/flowdiagnostics/AnisotropicEikonal.cpp
  opm/core/flowdiagnostics/DGBasis.cpp
  opm/core/flowdiagnostics/FlowDiagnostics.cpp
  opm/core/flowdiagnostics/TofDiscGalReorder.cpp
  opm/core/flowdiagnostics/TofReorder.cpp
  opm/core/linalg/LinearSolverFactory.cpp
  opm/core/linalg/LinearSolverInterface.cpp
  opm/core/linalg/LinearSolverIstl.cpp
  opm/core/linalg/LinearSolverPetsc.cpp
  opm/core/linalg/LinearSolverUmfpack.cpp
  opm/core/linalg/call_umfpack.c
  opm/core/linalg/sparse_sys.c
  opm/core/pressure/CompressibleTpfa.cpp
  opm/core/pressure/FlowBCManager.cpp
  opm/core/pressure/IncompTpfa.cpp
  opm/core/pressure/IncompTpfaSinglePhase.cpp
  opm/core/pressure/flow_bc.c
  opm/core/pressure/mimetic/mimetic.c
  opm/core/pressure/msmfem/dfs.c
  opm/core/pressure/msmfem/partition.c
  opm/core/pressure/tpfa/cfs_tpfa_residual.c
  opm/core/pressure/tpfa/ifs_tpfa.c
  opm/core/props/BlackoilPropertiesBasic.cpp
  opm/core/props/BlackoilPropertiesFromDeck.cpp
  opm/core/props/IncompPropertiesBasic.cpp
  opm/core/props/IncompPropertiesFromDeck.cpp
  opm/core/props/IncompPropertiesSinglePhase.cpp
  opm/core/props/pvt/PvtPropertiesBasic.cpp
  opm/core/props/pvt/PvtPropertiesIncompFromDeck.cpp
  opm/core/props/rock/RockBasic.cpp
  opm/core/props/rock/RockCompressibility.cpp
  opm/core/props/rock/RockFromDeck.cpp
  opm/core/props/satfunc/RelpermDiagnostics.cpp
  opm/core/props/satfunc/SaturationPropsBasic.cpp
  opm/core/props/satfunc/SaturationPropsFromDeck.cpp
  opm/core/simulator/BlackoilState.cpp
  opm/core/simulator/TwophaseState.cpp
  opm/core/simulator/SimulatorReport.cpp
  opm/core/transport/TransportSolverTwophaseInterface.cpp
  opm/core/transport/reorder/ReorderSolverInterface.cpp
  opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.cpp
  opm/core/transport/reorder/TransportSolverTwophaseReorder.cpp
  opm/core/transport/reorder/reordersequence.cpp
  opm/core/transport/reorder/tarjan.c
  opm/core/utility/Event.cpp
  opm/core/utility/miscUtilities.cpp
  opm/core/utility/miscUtilitiesBlackoil.cpp
  opm/core/utility/NullStream.cpp
  opm/core/wells/InjectionSpecification.cpp
  opm/core/wells/ProductionSpecification.cpp
  opm/core/wells/WellCollection.cpp
  opm/core/wells/WellsGroup.cpp
  opm/core/wells/WellsManager.cpp
  opm/core/wells/well_controls.c
  opm/core/wells/wells.c
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
  opm/simulators/ensureDirectoryExists.cpp
  opm/simulators/SimulatorCompressibleTwophase.cpp
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
  tests/test_multiphaseupwind.cpp
  tests/test_wellmodel.cpp
#  tests/test_thresholdpressure.cpp
  tests/test_wellswitchlogger.cpp
  tests/test_timer.cpp
  tests/test_invert.cpp
  tests/test_event.cpp
  tests/test_dgbasis.cpp
  tests/test_flowdiagnostics.cpp
  tests/test_parallelistlinformation.cpp
  tests/test_wells.cpp
  tests/test_linearsolver.cpp
  tests/test_parallel_linearsolver.cpp
  tests/test_satfunc.cpp
  tests/test_shadow.cpp
  tests/test_equil_legacy.cpp
  tests/test_blackoilstate.cpp
  tests/test_wellsmanager.cpp
  tests/test_wellcontrols.cpp
  tests/test_wellsgroup.cpp
  tests/test_wellcollection.cpp
  tests/test_pinchprocessor.cpp
  tests/test_anisotropiceikonal.cpp
  tests/test_stoppedwells.cpp
  tests/test_relpermdiagnostics.cpp
  tests/test_norne_pvt.cpp
  )

list (APPEND TEST_DATA_FILES
  tests/fluid.data
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
  tests/satfuncStandard.DATA
  tests/satfuncEPSBase.DATA
  tests/satfuncEPS_A.DATA
  tests/satfuncEPS_B.DATA
  tests/satfuncEPS_C.DATA
  tests/satfuncEPS_D.DATA
  tests/testBlackoilState1.DATA
  tests/testBlackoilState2.DATA
  tests/testPinch1.DATA
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
  examples/flow_legacy.cpp
  examples/flow_reorder.cpp
  examples/flow_sequential.cpp
  examples/flow.cpp
  examples/sim_2p_incomp_ad.cpp
  examples/sim_2p_comp_reorder.cpp
  examples/sim_simple.cpp
  examples/sim_poly2p_comp_reorder.cpp
  examples/sim_poly2p_incomp_reorder.cpp
  examples/wells_example.cpp
  examples/compute_eikonal_from_files.cpp
  examples/compute_initial_state.cpp
  examples/compute_tof_from_files.cpp
  examples/diagnose_relperm.cpp
  tutorials/sim_tutorial1.cpp
  tutorials/sim_tutorial2.cpp
  tutorials/sim_tutorial3.cpp
  tutorials/sim_tutorial4.cpp
  )

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
  examples/sim_2p_incomp.cpp
  examples/sim_2p_incomp_ad.cpp
  examples/sim_2p_comp_reorder.cpp
  examples/flow.cpp
  examples/flow_legacy.cpp
  examples/flow_reorder.cpp
  examples/flow_sequential.cpp
  examples/sim_poly2p_comp_reorder.cpp
  examples/sim_poly2p_incomp_reorder.cpp
  )

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/autodiff/AdditionalObjectDeleter.hpp
  opm/autodiff/AutoDiffBlock.hpp
  opm/autodiff/AutoDiffHelpers.hpp
  opm/autodiff/AutoDiffMatrix.hpp
  opm/autodiff/AutoDiff.hpp
  opm/autodiff/BlackoilDetails.hpp
  opm/autodiff/BlackoilLegacyDetails.hpp
  opm/autodiff/BlackoilModel.hpp
  opm/autodiff/BlackoilModelBase.hpp
  opm/autodiff/BlackoilModelBase_impl.hpp
  opm/autodiff/BlackoilModelEnums.hpp
  opm/autodiff/BlackoilModelParameters.hpp
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
  opm/autodiff/FlowMain.hpp
  opm/autodiff/FlowMainEbos.hpp
  opm/autodiff/FlowMainPolymer.hpp
  opm/autodiff/FlowMainSequential.hpp
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
  opm/autodiff/SimulatorIncompTwophaseAd.hpp
  opm/autodiff/SimulatorSequentialBlackoil.hpp
  opm/autodiff/TransportSolverTwophaseAd.hpp
  opm/autodiff/WellDensitySegmented.hpp
  opm/autodiff/WellStateFullyImplicitBlackoil.hpp
  opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp
  opm/autodiff/BlackoilOutputEbos.hpp
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
  opm/autodiff/ThreadHandle.hpp
  opm/core/flowdiagnostics/AnisotropicEikonal.hpp
  opm/core/flowdiagnostics/DGBasis.hpp
  opm/core/flowdiagnostics/FlowDiagnostics.hpp
  opm/core/flowdiagnostics/TofDiscGalReorder.hpp
  opm/core/flowdiagnostics/TofReorder.hpp
  opm/core/linalg/LinearSolverFactory.hpp
  opm/core/linalg/LinearSolverInterface.hpp
  opm/core/linalg/LinearSolverIstl.hpp
  opm/core/linalg/LinearSolverPetsc.hpp
  opm/core/linalg/LinearSolverUmfpack.hpp
  opm/core/linalg/ParallelIstlInformation.hpp
  opm/core/linalg/call_umfpack.h
  opm/core/linalg/sparse_sys.h
  opm/core/pressure/CompressibleTpfa.hpp
  opm/core/pressure/FlowBCManager.hpp
  opm/core/pressure/IncompTpfa.hpp
  opm/core/pressure/flow_bc.h
  opm/core/pressure/legacy_well.h
  opm/core/pressure/mimetic/mimetic.h
  opm/core/pressure/msmfem/dfs.h
  opm/core/pressure/msmfem/partition.h
  opm/core/pressure/tpfa/cfs_tpfa_residual.h
  opm/core/pressure/tpfa/compr_quant_general.h
  opm/core/pressure/tpfa/compr_source.h
  opm/core/pressure/tpfa/ifs_tpfa.h
  opm/core/props/BlackoilPhases.hpp
  opm/core/props/BlackoilPropertiesBasic.hpp
  opm/core/props/BlackoilPropertiesFromDeck.hpp
  opm/core/props/BlackoilPropertiesInterface.hpp
  opm/core/props/IncompPropertiesBasic.hpp
  opm/core/props/IncompPropertiesFromDeck.hpp
  opm/core/props/IncompPropertiesInterface.hpp
  opm/core/props/IncompPropertiesShadow.hpp
  opm/core/props/IncompPropertiesShadow_impl.hpp
  opm/core/props/IncompPropertiesSinglePhase.hpp
  opm/core/props/phaseUsageFromDeck.hpp
  opm/core/props/pvt/PvtPropertiesBasic.hpp
  opm/core/props/pvt/PvtPropertiesIncompFromDeck.hpp
  opm/core/props/pvt/ThermalGasPvtWrapper.hpp
  opm/core/props/pvt/ThermalOilPvtWrapper.hpp
  opm/core/props/pvt/ThermalWaterPvtWrapper.hpp
  opm/core/props/rock/RockBasic.hpp
  opm/core/props/rock/RockCompressibility.hpp
  opm/core/props/rock/RockFromDeck.hpp
  opm/core/props/satfunc/RelpermDiagnostics.hpp
  opm/core/props/satfunc/SaturationPropsBasic.hpp
  opm/core/props/satfunc/SaturationPropsFromDeck.hpp
  opm/core/props/satfunc/SaturationPropsInterface.hpp
  opm/core/props/satfunc/RelpermDiagnostics_impl.hpp
  opm/core/simulator/BlackoilState.hpp
  opm/core/simulator/BlackoilStateToFluidState.hpp
  opm/core/simulator/EquilibrationHelpers.hpp
  opm/core/simulator/ExplicitArraysFluidState.hpp
  opm/core/simulator/ExplicitArraysSatDerivativesFluidState.hpp
  opm/core/simulator/SimulatorReport.hpp
  opm/core/simulator/TwophaseState.hpp
  opm/core/simulator/WellState.hpp
  opm/core/simulator/initState.hpp
  opm/core/simulator/initStateEquil.hpp
  opm/core/simulator/initStateEquil_impl.hpp
  opm/core/simulator/initState_impl.hpp
  opm/core/transport/TransportSolverTwophaseInterface.hpp
  opm/core/transport/reorder/ReorderSolverInterface.hpp
  opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp
  opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp
  opm/core/transport/reorder/reordersequence.h
  opm/core/transport/reorder/tarjan.h
  opm/core/utility/DataMap.hpp
  opm/core/utility/Event.hpp
  opm/core/utility/Event_impl.hpp
  opm/core/utility/initHydroCarbonState.hpp
  opm/core/utility/miscUtilities.hpp
  opm/core/utility/miscUtilitiesBlackoil.hpp
  opm/core/utility/miscUtilities_impl.hpp
  opm/core/utility/NullStream.hpp
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
  opm/simulators/flow_ebos_blackoil.hpp
  opm/simulators/flow_ebos_gasoil.hpp
  opm/simulators/flow_ebos_oilwater.hpp
  opm/simulators/flow_ebos_polymer.hpp
  opm/simulators/flow_ebos_solvent.hpp
  opm/simulators/ensureDirectoryExists.hpp
  opm/simulators/ParallelFileMerger.hpp
  opm/simulators/SimulatorCompressibleTwophase.hpp
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
