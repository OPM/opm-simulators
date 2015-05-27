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

set( OPM_CMAKE_ROOT "" CACHE PATH "Root directory containing OPM related cmake modules")
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

if(NOT OPM_CMAKE_ROOT)
   find_package(opm-cmake QUIET)
endif()

if (opm-cmake_FOUND)
   include(OpmInit)
else()
   if (NOT OPM_CMAKE_ROOT AND SIBLING_SEARCH)
      set(OPM_CMAKE_ROOT ${PROJECT_SOURCE_DIR}/../opm-cmake)
   endif()
   if (OPM_CMAKE_ROOT)
      list( APPEND CMAKE_MODULE_PATH "${OPM_CMAKE_ROOT}/cmake/Modules")
      include (OpmInit OPTIONAL RESULT_VARIABLE OPM_INIT)
      set( OPM_MACROS_ROOT ${OPM_CMAKE_ROOT} )
   endif()

   if (NOT OPM_INIT)
      message( "" )
      message( " /-------------------------------------------------------------------------------\\")
      message( " |  Could not locate the opm build macros. The opm build macros                  |")
      message( " |  are in a separate repository - instructions to proceed:                      |")
      message( " |                                                                               |")
      message( " |    1. Clone the repository: git clone git@github.com:OPM/opm-cmake.git        |")
      message( " |                                                                               |")
      message( " |    2. Run cmake in the current project with -DOPM_CMAKE_ROOT=<path>/opm-cmake |")
      message( " |                                                                               |")
      message( " \\-------------------------------------------------------------------------------/")
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

add_custom_target(check-commits
									COMMAND ${CMAKE_COMMAND}
													-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
													-DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
													-P ${PROJECT_SOURCE_DIR}/cmake/Scripts/CheckCommits.cmake)
