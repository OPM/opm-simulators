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

  # these solvers are only compiled in if their dependency is found
  if (NOT dune-istl_FOUND)
    list (REMOVE_ITEM opm-simulators_SOURCES
      ${PROJECT_SOURCE_DIR}/opm/core/linalg/LinearSolverIstl.cpp
      )
  endif ()

  if (NOT PETSc_FOUND)
    list (REMOVE_ITEM opm-simulators_SOURCES
      ${PROJECT_SOURCE_DIR}/opm/core/linalg/LinearSolverPetsc.cpp
      )
  endif ()

  if ((NOT MPI_FOUND) OR (NOT DUNE_ISTL_FOUND))
    list (REMOVE_ITEM tests_SOURCES
      ${PROJECT_SOURCE_DIR}/tests/test_parallel_linearsolver.cpp
      ${PROJECT_SOURCE_DIR}/tests/test_parallelistlinformation.cpp
      )
  endif ()
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
