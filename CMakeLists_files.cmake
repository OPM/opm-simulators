# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
#	MAIN_SOURCE_FILES     List of compilation units which will be included in
#	                      the library. If it isn't on this list, it won't be
#	                      part of the library. Please try to keep it sorted to
#	                      maintain sanity.
#
#	TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
#	TEST_DATA_FILES       Files from the source three that should be made
#	                      available in the corresponding location in the build
#	                      tree in order to run tests there.
#
#	EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#	                      build, but which is not part of the library nor is
#	                      run as tests.
#
#	PUBLIC_HEADER_FILES   List of public header files that should be
#	                      distributed together with the library. The source
#	                      files can of course include other files than these;
#	                      you should only add to this list if the *user* of
#	                      the library needs it.

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	opm/core/grid/GridManager.cpp
	opm/core/grid/grid.c
	opm/core/grid/cart_grid.c
	opm/core/grid/cornerpoint_grid.c
	opm/core/grid/cpgpreprocess/facetopology.c
	opm/core/grid/cpgpreprocess/geometry.c
	opm/core/grid/cpgpreprocess/preprocess.c
	opm/core/grid/cpgpreprocess/uniquepoints.c
	opm/core/io/eclipse/EclipseGridInspector.cpp
	opm/core/io/eclipse/EclipseGridParser.cpp
	opm/core/io/eclipse/writeECLData.cpp
	opm/core/io/vag/vag.cpp
	opm/core/io/vtk/writeVtkData.cpp
	opm/core/linalg/LinearSolverAGMG.cpp
	opm/core/linalg/LinearSolverFactory.cpp
	opm/core/linalg/LinearSolverInterface.cpp
	opm/core/linalg/LinearSolverIstl.cpp
	opm/core/linalg/LinearSolverUmfpack.cpp
	opm/core/linalg/call_umfpack.c
	opm/core/linalg/sparse_sys.c
	opm/core/pressure/CompressibleTpfa.cpp
	opm/core/pressure/FlowBCManager.cpp
	opm/core/pressure/IncompTpfa.cpp
	opm/core/pressure/cfsh.c
	opm/core/pressure/flow_bc.c
	opm/core/pressure/fsh.c
	opm/core/pressure/fsh_common_impl.c
	opm/core/pressure/ifsh.c
	opm/core/pressure/mimetic/hybsys.c
	opm/core/pressure/mimetic/hybsys_global.c
	opm/core/pressure/mimetic/mimetic.c
	opm/core/pressure/msmfem/coarse_conn.c
	opm/core/pressure/msmfem/coarse_sys.c
	opm/core/pressure/msmfem/dfs.c
	opm/core/pressure/msmfem/hash_set.c
	opm/core/pressure/msmfem/ifsh_ms.c
	opm/core/pressure/msmfem/partition.c
	opm/core/pressure/tpfa/cfs_tpfa.c
	opm/core/pressure/tpfa/cfs_tpfa_residual.c
	opm/core/pressure/tpfa/compr_bc.c
	opm/core/pressure/tpfa/compr_quant.c
	opm/core/pressure/tpfa/compr_quant_general.c
	opm/core/pressure/tpfa/compr_source.c
	opm/core/pressure/tpfa/ifs_tpfa.c
	opm/core/pressure/tpfa/trans_tpfa.c
	opm/core/pressure/legacy_well.c
	opm/core/props/BlackoilPropertiesBasic.cpp
	opm/core/props/BlackoilPropertiesFromDeck.cpp
	opm/core/props/IncompPropertiesBasic.cpp
	opm/core/props/IncompPropertiesFromDeck.cpp
	opm/core/props/pvt/BlackoilPvtProperties.cpp
	opm/core/props/pvt/PvtPropertiesBasic.cpp
	opm/core/props/pvt/PvtPropertiesIncompFromDeck.cpp
	opm/core/props/pvt/SinglePvtDead.cpp
	opm/core/props/pvt/SinglePvtDeadSpline.cpp
	opm/core/props/pvt/SinglePvtInterface.cpp
	opm/core/props/pvt/SinglePvtLiveGas.cpp
	opm/core/props/pvt/SinglePvtLiveOil.cpp
	opm/core/props/rock/RockBasic.cpp
	opm/core/props/rock/RockCompressibility.cpp
	opm/core/props/rock/RockFromDeck.cpp
	opm/core/props/satfunc/SatFuncGwseg.cpp
	opm/core/props/satfunc/SatFuncSimple.cpp
	opm/core/props/satfunc/SatFuncStone2.cpp
	opm/core/props/satfunc/SaturationPropsBasic.cpp
	opm/core/props/satfunc/SaturationPropsFromDeck.cpp
	opm/core/simulator/SimulatorCompressibleTwophase.cpp
	opm/core/simulator/SimulatorIncompTwophase.cpp
	opm/core/simulator/SimulatorReport.cpp
	opm/core/simulator/SimulatorTimer.cpp
	opm/core/tof/DGBasis.cpp
	opm/core/tof/TofReorder.cpp
	opm/core/tof/TofDiscGalReorder.cpp
	opm/core/transport/TransportSolverTwophaseInterface.cpp
	opm/core/transport/implicit/TransportSolverTwophaseImplicit.cpp
	opm/core/transport/implicit/transport_source.c
	opm/core/transport/minimal/spu_explicit.c
	opm/core/transport/minimal/spu_implicit.c
	opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.cpp
	opm/core/transport/reorder/ReorderSolverInterface.cpp
	opm/core/transport/reorder/TransportSolverTwophaseReorder.cpp
	opm/core/transport/reorder/reordersequence.cpp
	opm/core/transport/reorder/tarjan.c
	opm/core/utility/MonotCubicInterpolator.cpp
	opm/core/utility/StopWatch.cpp
	opm/core/utility/VelocityInterpolation.cpp
	opm/core/utility/WachspressCoord.cpp
	opm/core/utility/miscUtilities.cpp
	opm/core/utility/miscUtilitiesBlackoil.cpp
	opm/core/utility/parameters/Parameter.cpp
	opm/core/utility/parameters/ParameterGroup.cpp
	opm/core/utility/parameters/ParameterTools.cpp
	opm/core/utility/parameters/ParameterXML.cpp
	opm/core/utility/parameters/tinyxml/tinystr.cpp
	opm/core/utility/parameters/tinyxml/tinyxml.cpp
	opm/core/utility/parameters/tinyxml/tinyxmlerror.cpp
	opm/core/utility/parameters/tinyxml/tinyxmlparser.cpp
	opm/core/utility/parameters/tinyxml/xmltest.cpp
	opm/core/wells/InjectionSpecification.cpp
	opm/core/wells/ProductionSpecification.cpp
	opm/core/wells/WellCollection.cpp
	opm/core/wells/WellsGroup.cpp
	opm/core/wells/WellsManager.cpp
	opm/core/wells/wells.c
	)

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
	tests/test_dgbasis.cpp
	tests/test_nonuniformtablelinear.cpp
	tests/test_sparsevector.cpp
	tests/test_sparsetable.cpp
	tests/test_velocityinterpolation.cpp
	tests/test_quadratures.cpp
	tests/test_uniformtablelinear.cpp
	tests/test_wells.cpp
	tests/test_wachspresscoord.cpp
	tests/test_column_extract.cpp
	tests/test_geom2d.cpp
	tests/test_param.cpp
	)

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
	tests/extratestdata.xml
	tests/testdata.xml
	)

# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/compute_tof.cpp
	examples/compute_tof_from_files.cpp
	examples/import_rewrite.cpp
	examples/sim_2p_comp_reorder.cpp
	examples/sim_2p_incomp.cpp
	examples/wells_example.cpp
	tutorials/tutorial1.cpp
	tutorials/tutorial2.cpp
	tutorials/tutorial3.cpp
	tutorials/tutorial4.cpp
	)

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
	examples/sim_2p_comp_reorder.cpp
	examples/sim_2p_incomp.cpp
	)

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
	opm/core/doxygen_main.hpp
	opm/core/grid.h
	opm/core/grid/CellQuadrature.hpp
	opm/core/grid/ColumnExtract.hpp
	opm/core/grid/FaceQuadrature.hpp
	opm/core/grid/GridManager.hpp
	opm/core/grid/cart_grid.h
	opm/core/grid/cornerpoint_grid.h
	opm/core/grid/cpgpreprocess/facetopology.h
	opm/core/grid/cpgpreprocess/geometry.h
	opm/core/grid/cpgpreprocess/grdecl.h
	opm/core/grid/cpgpreprocess/preprocess.h
	opm/core/grid/cpgpreprocess/uniquepoints.h
	opm/core/io/eclipse/CornerpointChopper.hpp
	opm/core/io/eclipse/EclipseGridInspector.hpp
	opm/core/io/eclipse/EclipseGridParser.hpp
	opm/core/io/eclipse/EclipseGridParserHelpers.hpp
	opm/core/io/eclipse/EclipseUnits.hpp
	opm/core/io/eclipse/SpecialEclipseFields.hpp
	opm/core/io/eclipse/writeECLData.hpp
	opm/core/io/vag/vag.hpp
	opm/core/io/vtk/writeVtkData.hpp
	opm/core/linalg/LinearSolverAGMG.hpp
	opm/core/linalg/LinearSolverFactory.hpp
	opm/core/linalg/LinearSolverInterface.hpp
	opm/core/linalg/LinearSolverIstl.hpp
	opm/core/linalg/LinearSolverUmfpack.hpp
	opm/core/linalg/blas_lapack.h
	opm/core/linalg/call_umfpack.h
	opm/core/linalg/sparse_sys.h
	opm/core/wells.h
	opm/core/pressure/CompressibleTpfa.hpp
	opm/core/pressure/FlowBCManager.hpp
	opm/core/pressure/IncompTpfa.hpp
	opm/core/pressure/flow_bc.h
	opm/core/pressure/fsh.h
	opm/core/pressure/fsh_common_impl.h
	opm/core/pressure/legacy_well.h
	opm/core/pressure/mimetic/hybsys.h
	opm/core/pressure/mimetic/hybsys_global.h
	opm/core/pressure/mimetic/mimetic.h
	opm/core/pressure/msmfem/coarse_conn.h
	opm/core/pressure/msmfem/coarse_sys.h
	opm/core/pressure/msmfem/dfs.h
	opm/core/pressure/msmfem/hash_set.h
	opm/core/pressure/msmfem/ifsh_ms.h
	opm/core/pressure/msmfem/partition.h
	opm/core/pressure/tpfa/cfs_tpfa.h
	opm/core/pressure/tpfa/cfs_tpfa_residual.h
	opm/core/pressure/tpfa/compr_bc.h
	opm/core/pressure/tpfa/compr_quant.h
	opm/core/pressure/tpfa/compr_quant_general.h
	opm/core/pressure/tpfa/compr_source.h
	opm/core/pressure/tpfa/ifs_tpfa.h
	opm/core/pressure/tpfa/trans_tpfa.h
	opm/core/props/BlackoilPhases.hpp
	opm/core/props/BlackoilPropertiesBasic.hpp
	opm/core/props/BlackoilPropertiesFromDeck.hpp
	opm/core/props/BlackoilPropertiesInterface.hpp
	opm/core/props/IncompPropertiesBasic.hpp
	opm/core/props/IncompPropertiesFromDeck.hpp
	opm/core/props/IncompPropertiesInterface.hpp
	opm/core/props/phaseUsageFromDeck.hpp
	opm/core/props/pvt/BlackoilPvtProperties.hpp
	opm/core/props/pvt/PvtPropertiesBasic.hpp
	opm/core/props/pvt/PvtPropertiesIncompFromDeck.hpp
	opm/core/props/pvt/SinglePvtConstCompr.hpp
	opm/core/props/pvt/SinglePvtDead.hpp
	opm/core/props/pvt/SinglePvtDeadSpline.hpp
	opm/core/props/pvt/SinglePvtInterface.hpp
	opm/core/props/pvt/SinglePvtLiveGas.hpp
	opm/core/props/pvt/SinglePvtLiveOil.hpp
	opm/core/props/rock/RockBasic.hpp
	opm/core/props/rock/RockCompressibility.hpp
	opm/core/props/rock/RockFromDeck.hpp
	opm/core/props/satfunc/SatFuncGwseg.hpp
	opm/core/props/satfunc/SatFuncSimple.hpp
	opm/core/props/satfunc/SatFuncStone2.hpp
	opm/core/props/satfunc/SaturationPropsBasic.hpp
	opm/core/props/satfunc/SaturationPropsFromDeck.hpp
	opm/core/props/satfunc/SaturationPropsFromDeck_impl.hpp
	opm/core/props/satfunc/SaturationPropsInterface.hpp
	opm/core/simulator/BlackoilState.hpp
	opm/core/simulator/SimulatorCompressibleTwophase.hpp
	opm/core/simulator/SimulatorIncompTwophase.hpp
	opm/core/simulator/SimulatorReport.hpp
	opm/core/simulator/SimulatorTimer.hpp
	opm/core/simulator/TwophaseState.hpp
	opm/core/simulator/WellState.hpp
	opm/core/simulator/initState.hpp
	opm/core/simulator/initState_impl.hpp
	opm/core/tof/DGBasis.hpp
	opm/core/tof/TofReorder.hpp
	opm/core/tof/TofDiscGalReorder.hpp
	opm/core/transport/TransportSolverTwophaseInterface.hpp
	opm/core/transport/implicit/CSRMatrixBlockAssembler.hpp
	opm/core/transport/implicit/CSRMatrixUmfpackSolver.hpp
	opm/core/transport/implicit/GravityColumnSolver.hpp
	opm/core/transport/implicit/GravityColumnSolver_impl.hpp
	opm/core/transport/implicit/ImplicitAssembly.hpp
	opm/core/transport/implicit/ImplicitTransport.hpp
	opm/core/transport/implicit/JacobianSystem.hpp
	opm/core/transport/implicit/NormSupport.hpp
	opm/core/transport/implicit/SinglePointUpwindTwoPhase.hpp
	opm/core/transport/implicit/TransportSolverTwophaseImplicit.hpp
	opm/core/transport/implicit/SimpleFluid2pWrappingProps.hpp
	opm/core/transport/implicit/SimpleFluid2pWrappingProps_impl.hpp
	opm/core/transport/implicit/transport_source.h
	opm/core/transport/minimal/spu_explicit.h
	opm/core/transport/minimal/spu_implicit.h
	opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp
	opm/core/transport/reorder/ReorderSolverInterface.hpp
	opm/core/transport/reorder/TransportSolverTwophaseReorder.hpp
	opm/core/transport/reorder/reordersequence.h
	opm/core/transport/reorder/tarjan.h
	opm/core/utility/Average.hpp
	opm/core/utility/DataMap.hpp
	opm/core/utility/ErrorMacros.hpp
	opm/core/utility/Factory.hpp
	opm/core/utility/MonotCubicInterpolator.hpp
	opm/core/utility/NonuniformTableLinear.hpp
	opm/core/utility/RootFinders.hpp
	opm/core/utility/SparseTable.hpp
	opm/core/utility/SparseVector.hpp
	opm/core/utility/StopWatch.hpp
	opm/core/utility/UniformTableLinear.hpp
	opm/core/utility/Units.hpp
	opm/core/utility/VelocityInterpolation.hpp
	opm/core/utility/WachspressCoord.hpp
	opm/core/utility/buildUniformMonotoneTable.hpp
	opm/core/utility/have_boost_redef.hpp
	opm/core/utility/linearInterpolation.hpp
	opm/core/utility/miscUtilities.hpp
	opm/core/utility/miscUtilitiesBlackoil.hpp
	opm/core/utility/parameters/Parameter.hpp
	opm/core/utility/parameters/ParameterGroup.hpp
	opm/core/utility/parameters/ParameterGroup_impl.hpp
	opm/core/utility/parameters/ParameterMapItem.hpp
	opm/core/utility/parameters/ParameterRequirement.hpp
	opm/core/utility/parameters/ParameterStrings.hpp
	opm/core/utility/parameters/ParameterTools.hpp
	opm/core/utility/parameters/ParameterXML.hpp
	opm/core/utility/parameters/tinyxml/tinystr.h
	opm/core/utility/parameters/tinyxml/tinyxml.h
	opm/core/wells/InjectionSpecification.hpp
	opm/core/wells/ProductionSpecification.hpp
	opm/core/wells/WellCollection.hpp
	opm/core/wells/WellsGroup.hpp
	opm/core/wells/WellsManager.hpp
	)
