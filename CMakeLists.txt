# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

cmake_minimum_required(VERSION 3.0)

# specify the option that allows to select the build system. if
# dunecontrol is used, we default to the dune build system. if not, we
# use the OPM home-brewn one
set(USE_DUNE_BUILDSYSTEM_DEFAULT OFF)
if (DEFINED ENV{DUNE_CONTROL_PATH})
  set(USE_DUNE_BUILDSYSTEM_DEFAULT ON)
endif()
option(USE_DUNE_BUILDSYSTEM "Use the Dune build system?" ${USE_DUNE_BUILDSYSTEM_DEFAULT})

if (USE_DUNE_BUILDSYSTEM)

######
# use DUNE's build system
######

# set up project and specify the minimum cmake version
project("opm-simulators" C CXX)

# Sibling build
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

# guess the sibling dir
get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
foreach(MOD dune-common dune-geometry dune-grid dune-localfunctions dune-istl opm-common libecl opm-parser opm-grid opm-data opm-core opm-output opm-material ewoms)
	# do not overwrite the location of a module if it was explicitly
	# specified by the user
	if(${MOD}_DIR)
		continue()
	endif()

	# Try various possible locations for the build directory of the dependency
	foreach(BUILD_DIR "${_leaf_dir_name}" "build-cmake" "build" "../${MOD}-build" "../build-${MOD}" "../build" ".")
		if(EXISTS "${_modules_dir}/${MOD}/${BUILD_DIR}/${MOD}-config.cmake")
			set(${MOD}_DIR ${_modules_dir}/${MOD}/${BUILD_DIR} PARENT_SCOPE)
			break()
		endif()
	endforeach()
endforeach()

if(dune-common_DIR AND NOT IS_DIRECTORY ${dune-common_DIR})
  message(WARNING "Value ${dune-common_DIR} passed to variable"
    " dune-common_DIR is not a directory")
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()


# Set CMP0053 (how to handle escape sequences in strings) to the new
# behavior to avoid a pretty annoying cmake warning if a library is
# defined in the toplevel CMakeLists.txt. This should probably be
# considered to be a bug in the dune build system. Note that the old
# behaviour will most likely also work fine, but the result of setting
# this policy to NEW is most likely what's intended.
if (POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif()

# find the build system (i.e., dune-common) and set cmake's module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# include the OPM cmake macros
include(OpmMacros)

# find the packages needed to compile the module
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

find_package(PlainOpmData)
find_package(OpmOutputUtils)

# recursively mark all header files beneath the "opm" directory for
# installation.
opm_recusive_export_all_headers("opm")

# the cmake modules get a special treatment
opm_export_cmake_modules()

# we want all features detected by the build system to be enabled,
# thank you!
dune_enable_all_packages()

# add source files from CMakeLists_files.cmake to library and create executables
opm_add_library_and_executables("opmsimulators")

# download Eigen if user doesn't have the correct version
if (NOT EIGEN3_FOUND)
        message (STATUS "Downloading Eigen3")
        include (ExternalProject)
        externalProject_Add (Eigen3
                GIT_REPOSITORY git://github.com/OPM/eigen3
                UPDATE_COMMAND git checkout 9e788db99d73f3199ade74bfda8d9f73fdfbbe4c
                CMAKE_ARGS -Wno-dev -DEIGEN_TEST_NO_OPENGL=1 -DEIGEN_BUILD_PKGCONFIG=0 -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/eigen3-installed
                )

        include_directories (${CMAKE_BINARY_DIR}/eigen3-installed/include/eigen3)
        add_dependencies (opmsimulators Eigen3)
endif (NOT EIGEN3_FOUND)

if(opm-data_FOUND)
  include (${CMAKE_CURRENT_SOURCE_DIR}/compareECLFiles.cmake)
endif()

opm_recusive_copy_testdata("tests/*.data" "tests/*.DATA" "tests/VFP*")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

else()

######
# use the traditional OPM build system
######

# Mandatory call to project
project(opm-simulators CXX)

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental

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
      set(opm-common_DIR "${_parent_full_dir}/opm-common}")
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
  if (NOT SuiteSparse_FOUND)
    list (REMOVE_ITEM examples_SOURCES
      ${PROJECT_SOURCE_DIR}/tutorials/tutorial2.cpp
      ${PROJECT_SOURCE_DIR}/tutorials/tutorial3.cpp
      ${PROJECT_SOURCE_DIR}/tutorials/tutorial4.cpp
      )
  endif()
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
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
    CMAKE_ARGS -Wno-dev -DEIGEN_TEST_NO_OPENGL=1 -DEIGEN_BUILD_PKGCONFIG=0 -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/eigen3-installed
    )

  include_directories (${CMAKE_BINARY_DIR}/eigen3-installed/include/eigen3)
  add_dependencies (opmsimulators Eigen3)
endif (NOT EIGEN3_FOUND)



if (HAVE_OPM_DATA)
    include (${CMAKE_CURRENT_SOURCE_DIR}/compareECLFiles.cmake)
endif()

endif()
