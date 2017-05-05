# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

cmake_minimum_required (VERSION 2.8)

set( OPM_COMMON_ROOT "" CACHE PATH "Root directory containing OPM related cmake modules")
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental

if(NOT OPM_COMMON_ROOT)
   find_package(opm-common QUIET)
endif()

if (opm-common_FOUND)
   include(OpmInit)
else()
   unset(opm-common_FOUND)

   if (NOT OPM_COMMON_ROOT AND SIBLING_SEARCH)
      set(OPM_COMMON_ROOT ${PROJECT_SOURCE_DIR}/../opm-common)
   endif()
   if (OPM_COMMON_ROOT)
      list( APPEND CMAKE_MODULE_PATH "${OPM_COMMON_ROOT}/cmake/Modules")
      include (OpmInit OPTIONAL RESULT_VARIABLE OPM_INIT)
      set( OPM_MACROS_ROOT ${OPM_COMMON_ROOT} )
   endif()

   if (NOT OPM_INIT)
      message( "" )
      message( " /---------------------------------------------------------------------------------\\")
      message( " |  Could not locate the opm build macros. The opm build macros                    |")
      message( " |  are in a separate repository - instructions to proceed:                        |")
      message( " |                                                                                 |")
      message( " |    1. Clone the repository: git clone git@github.com:OPM/opm-common.git         |")
      message( " |                                                                                 |")
      message( " |    2. Run cmake in the current project with -DOPM_COMMON_ROOT=<path>/opm-common |")
      message( " |                                                                                 |")
      message( " \\---------------------------------------------------------------------------------/")
      message( "" )
      message( FATAL_ERROR "Could not find OPM Macros")
   endif()

endif()


# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# Look for the opm-data repository; if found the variable
# HAVE_OPM_DATA will be set to true.
include( Findopm-data )

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include ("${project}-prereqs")

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
	opm_need_version_of ("dune-common")
	opm_need_version_of ("dune-istl")
	opm_need_version_of ("ewoms")
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
	if(OPM_GRID_FOUND OR opm-grid_FOUND)
		list (APPEND examples_SOURCES
			${PROJECT_SOURCE_DIR}/examples/flow_mpi.cpp
			${PROJECT_SOURCE_DIR}/examples/flow_multisegment_mpi.cpp
			)
		list (APPEND examples_SOURCES_DIST
			${PROJECT_SOURCE_DIR}/examples/flow_mpi.cpp
			)
	endif()
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
  if(RUN_INTERNAL_TESTS)
    if(ERT_PYTHON_PATH)
      include (ExternalProject)
      # Pull down the opm-data repository
      externalproject_add(opm-data
                          PREFIX tests/opm-data
                          GIT_REPOSITORY git@github.com:OPM/opm-data
                          GIT_TAG origin/master
                          CONFIGURE_COMMAND ""
                          BUILD_COMMAND ""
                          INSTALL_COMMAND "")
      set_target_properties(opm-data PROPERTIES EXCLUDE_FROM_ALL 1)

      set(INPUT_DATA_PATH ${CMAKE_BINARY_DIR}/tests/opm-data/src/opm-data)

      # Add SPE1 test
      opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/runtest.sh "")
      opm_add_test(SPE1 NO_COMPILE
                   EXE_NAME ${CMAKE_BINARY_DIR}/bin/flow
                   DRIVER_ARGS ${PYTHON_EXECUTABLE}
                               ${ERT_PYTHON_PATH}
                               ${INPUT_DATA_PATH}/output_comparator/src/compare_eclipse.py
                               ${INPUT_DATA_PATH}/spe1/eclipse-simulation/
                               ${INPUT_DATA_PATH}/spe1/opm-simulation/
                               SPE1DECK 0.0006 0.004
                   TEST_ARGS ${INPUT_DATA_PATH}/spe1 spe1deck.xml)

      add_dependencies(check opm-data flow)
    else()
      message(WARNING "Internal tests cannot be enabled if ERT is built without Python support.")
    endif()
  endif()
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

# download Eigen if user doesn't have the correct version
if (NOT EIGEN3_FOUND)
	message (STATUS "Downloading Eigen3")
	include (ExternalProject)
	externalProject_Add (Eigen3
		GIT_REPOSITORY git://github.com/OPM/eigen3
		UPDATE_COMMAND git checkout 9e788db99d73f3199ade74bfda8d9f73fdfbbe4c
		CMAKE_ARGS -DEIGEN_TEST_NO_OPENGL=1 -DEIGEN_BUILD_PKGCONFIG=0 -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/eigen3-installed
		)

	include_directories (${CMAKE_BINARY_DIR}/eigen3-installed/include/eigen3)
	add_dependencies (opmsimulators Eigen3)
endif (NOT EIGEN3_FOUND)



if (HAVE_OPM_DATA)
    include (${CMAKE_CURRENT_SOURCE_DIR}/compareECLFiles.cmake)
endif()

# create a symbolic link from flow to flow_legacy
ADD_CUSTOM_TARGET(flow ALL
                  DEPENDS flow_ebos
                  COMMAND ${CMAKE_COMMAND} -E create_symlink "flow_ebos" "${CMAKE_BINARY_DIR}/bin/flow")
install(FILES ${CMAKE_BINARY_DIR}/bin/flow DESTINATION bin)
