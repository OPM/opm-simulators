# -*- mode: cmake; tab-width: 2; indent-tabs-mode: nil; truncate-lines: t; compile-command: "cmake -Wdev" -*-
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
project(opm-material CXX)

cmake_minimum_required (VERSION 2.8)

set( OPM_COMMON_ROOT "" CACHE PATH "Root directory containing OPM related cmake modules")
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

if(NOT OPM_COMMON_ROOT)
  if(SIBLING_SEARCH AND NOT opm-common_DIR)
    # guess the sibling dir
    get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
    get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
    get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
    if(_parent_dir_name STREQUAL CMAKE_PROJECT_NAME)
      #looks like <module-name>/<build-dir> is used
      get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
      if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
        set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
      endif()
    elseif(_leaf_dir_name STREQUAL CMAKE_PROJECT_NAME)
      # we assume <build-dir>/<module-name>
      if(IS_DIRECTORY ${_parent_full_dir}/opm-common)
        set(opm-common_DIR "${_parent_full_dir}/opm-common}")
      endif()
    endif()
  endif()
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

# not the same location as most of the other projects; this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include (${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
opm_need_version_of ("dune-common")
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

opm_add_test(test_eclblackoilpvt CONDITION OPM_PARSER_FOUND)
opm_add_test(test_eclmateriallawmanager CONDITION OPM_PARSER_FOUND)
opm_add_test(test_fluidmatrixinteractions)
opm_add_test(test_pengrobinson)
opm_add_test(test_densead)
opm_add_test(test_ncpflash)
opm_add_test(test_spline)
opm_add_test(test_tabulation)
opm_add_test(test_2dtables)
opm_add_test(test_components)
opm_add_test(test_fluidsystems)
opm_add_test(test_immiscibleflash)
