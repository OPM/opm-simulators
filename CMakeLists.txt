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

# Mandatory call to project
project(opm-simulators C CXX)

cmake_minimum_required (VERSION 2.8)

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental
option(BUILD_FLOW "Build the production oriented flow simulator?" ON)
option(BUILD_EBOS "Build the research oriented ebos simulator?" ON)
option(BUILD_EBOS_EXTENSIONS "Build the variants for various extensions of ebos by default?" OFF)

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)


# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# Look for the opm-tests repository; if found the variable
# HAVE_OPM_TESTS will be set to true.
include(Findopm-tests)

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
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

if (HAVE_OPM_TESTS)
    include (${CMAKE_CURRENT_SOURCE_DIR}/compareECLFiles.cmake)
endif()

opm_set_test_driver(${CMAKE_CURRENT_SOURCE_DIR}/tests/run-parallel-unitTest.sh "")

opm_add_test(test_gatherconvergencereport
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_gatherconvergencereport.cpp
  CONDITION
    MPI_FOUND
  DRIVER_ARGS
    5 ${CMAKE_BINARY_DIR}
)

opm_add_test(test_gatherdeferredlogger
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_gatherdeferredlogger.cpp
  CONDITION
    MPI_FOUND
  DRIVER_ARGS
    5 ${CMAKE_BINARY_DIR}
)

include(OpmBashCompletion)

if (NOT BUILD_FLOW)
  set(FLOW_DEFAULT_ENABLE_IF "FALSE")
else()
  set(FLOW_DEFAULT_ENABLE_IF "TRUE")
endif()

# the production oriented general-purpose ECL simulator
opm_add_test(flow
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEFAULT_ENABLE_IF ${FLOW_DEFAULT_ENABLE_IF}
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators"
  SOURCES
  flow/flow.cpp
  flow/flow_ebos_blackoil.cpp
  flow/flow_ebos_gasoil.cpp
  flow/flow_ebos_oilwater.cpp
  flow/flow_ebos_polymer.cpp
  flow/flow_ebos_solvent.cpp
  flow/flow_ebos_energy.cpp
  flow/flow_ebos_oilwater_polymer.cpp
  flow/flow_ebos_oilwater_polymer_injectivity.cpp)

if (BUILD_FLOW)
  install(TARGETS flow DESTINATION bin)
  opm_add_bash_completion(flow)

  add_test(NAME flow__version
    COMMAND flow --version)
  set_tests_properties(flow__version PROPERTIES
    PASS_REGULAR_EXPRESSION "${${project}_LABEL}")
endif()

if (NOT BUILD_EBOS)
  set(EBOS_DEFAULT_ENABLE_IF "FALSE")
else()
  set(EBOS_DEFAULT_ENABLE_IF "TRUE")
endif()
  
# the research oriented general-purpose ECL simulator ("ebos" == &ecl
# &black-&oil &simulator)
opm_add_test(ebos
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_DEFAULT_ENABLE_IF}
  ALWAYS_ENABLE
  EXE_NAME ebos
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators"
  SOURCES ebos/ebos.cc)

if (BUILD_EBOS)
  install(TARGETS ebos DESTINATION bin)
  opm_add_bash_completion(ebos)
endif()

if (NOT BUILD_EBOS_EXTENSIONS)
  set(EBOS_EXTENSIONS_DEFAULT_ENABLE_IF "FALSE")
else()
  set(EBOS_EXTENSIONS_DEFAULT_ENABLE_IF "TRUE")
endif()

opm_add_test(ebos_solvent
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/ebos_solvent.cc
  EXE_NAME ebos_solvent
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators")

opm_add_test(ebos_polymer
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/ebos_polymer.cc
  EXE_NAME ebos_polymer
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators")

opm_add_test(ebos_gasoil
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/ebos_gasoil.cc
  EXE_NAME ebos_gasoil
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators")

opm_add_test(ebos_oilwater
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/ebos_oilwater.cc
  EXE_NAME ebos_oilwater
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators")

opm_add_test(ebos_thermal
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/ebos_thermal.cc
  EXE_NAME ebos_thermal
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators")

opm_add_bash_completion(ebos_solvent)
opm_add_bash_completion(ebos_polymer)
opm_add_bash_completion(ebos_gasoil)
opm_add_bash_completion(ebos_oilwater)
opm_add_bash_completion(ebos_thermal)

if (BUILD_EBOS_EXTENSIONS)
  install(TARGETS ebos_solvent DESTINATION bin)
  install(TARGETS ebos_polymer DESTINATION bin)
  install(TARGETS ebos_gasoil DESTINATION bin)
  install(TARGETS ebos_oilwater DESTINATION bin)
  install(TARGETS ebos_thermal DESTINATION bin)
endif()

