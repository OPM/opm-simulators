# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# key information about the library
set (project "opm-core")
set (${project}_NAME "${project}")
set (${project}_DESCRIPTION "Open Porous Media Initiative Core Library")
set (${project}_DIR "opm")
set (${project}_VERSION_MAJOR 1)
set (${project}_VERSION_MINOR 0)
set (doxy_dir "Documentation")

# defines that must be present in config.h for our headers
set (${project}_CONFIG_VAR
	HAVE_AGMG
	HAVE_DUNE_ISTL
	HAVE_DYNAMIC_BOOST_TEST
	HAVE_ERT
	HAVE_SUITESPARSE_UMFPACK_H
	HAVE_NULLPTR
	HAVE_STATIC_ASSERT
	)

# dependencies
set (${project}_DEPS
	# compile with C99 support if available
	"C99"
	# compile with C++0x/11 support if available
	"CXX11Features"
	# various runtime library enhancements
	"Boost 1.39.0
		COMPONENTS date_time filesystem system unit_test_framework REQUIRED"
	# matrix library
	"BLAS REQUIRED"
	"LAPACK REQUIRED"
	# Tim Davis' SuiteSparse archive
	"SuiteSparse COMPONENTS umfpack"
	# solver
	"SuperLU"
	# xml processing (for config parsing)
	"TinyXML"
	# Ensembles-based Reservoir Tools (ERT)
	"ERT"
	# DUNE dependency
	"dune-common"
	"dune-istl"
	)

# additional search modules
set (${project}_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
list (APPEND CMAKE_MODULE_PATH ${${project}_MODULE_DIR})

# include special
if (CMAKE_VERSION VERSION_LESS "2.8.5")
	message (STATUS "Enabling compatibility modules for CMake 2.8.5")
	list (APPEND CMAKE_MODULE_PATH "${${project}_MODULE_DIR}/compat-2.8.5")	
endif (CMAKE_VERSION VERSION_LESS "2.8.5")	

if (CMAKE_VERSION VERSION_LESS "2.8.7")
	message (STATUS "Enabling compatibility modules for CMake 2.8.7")
	list (APPEND CMAKE_MODULE_PATH "${${project}_MODULE_DIR}/compat-2.8.7")
endif (CMAKE_VERSION VERSION_LESS "2.8.7")

# don't write default flags into the cache, preserve that for user set values
include (AddOptions)
no_default_options ()

# C++ project
cmake_minimum_required (VERSION 2.8)
project (${${project}_NAME})
enable_language (C)
enable_language (CXX)

# print system information to better pinpoint issues from log alone
include (UseSystemInfo)
system_info ()

# very early try to print repo id (to pinpoint version if something goes wrong)
include (UseVCSInfo)
vcs_info ()

# print toolchain information to identify compilers with potential bugs
include (UseCompVer)
compiler_info ()

# default settings: build static debug library
include (OpmDefaults)
opm_defaults (${project})
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# use tricks to do faster builds
include (UseFastBuilds)

# precompiled headers
include (UsePrecompHeaders)

# optimize full if we're not doing a debug build
include (UseOptimization)

# turn on all warnings; this must be done before adding any
# dependencies, in case they alter the list of warnings
include (UseWarnings)

# parallel computing must be explicitly enabled
option (USE_MPI "Use Message Passing Interface for parallel computing" OFF)
if (NOT USE_MPI)
	set (CMAKE_DISABLE_FIND_PACKAGE_MPI TRUE)
endif (NOT USE_MPI)

### --- begin opm-core specific --- ###
include (UseOpenMP)
find_openmp (${project})
### --- end opm-core specific --- ###

# macro to set standard variables (INCLUDE_DIRS, LIBRARIES etc.)
include (OpmFind)
find_and_append_package_list_to (${project} ${${project}_DEPS})

# remove the dependency on the testing framework from the main library;
# it is not possible to query for Boost twice with different components.
list (REMOVE_ITEM ${project}_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# don't import more libraries than we need to
include (UseOnlyNeeded)

# put debug information into every executable
include (UseDebugSymbols)

# detect if Boost is in a shared library
include (UseDynamicBoost)

# needed for Debian installation scheme
include (UseMultiArch)

# this module contains code to figure out which files is where
include (OpmFiles)
opm_auto_dirs ()

# put libraries in lib/
opm_out_dirs ()

# identify the compilation units in the library; sources in opm/,
# tests files in tests/, examples in tutorials/ and examples/
opm_sources (${project})

### --- begin AGMG specific --- ###
# Algebraic Multigrid must be compiled together with our program;
# if it is not available, then remove our corresponding component
find_package (AGMG)
if (AGMG_FOUND)
	list (APPEND ${project}_SOURCES ${AGMG_SOURCES})
endif (AGMG_FOUND)
### --- end AGMG specific --- ###

### --- begin opm-core specific --- ###
# these solvers are only compiled in if their dependency is found
if (NOT AGMG_FOUND)
	list (REMOVE_ITEM opm-core_SOURCES
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/LinearSolverAGMG.cpp
		)
endif (NOT AGMG_FOUND)
if (NOT dune-istl_FOUND)
	list (REMOVE_ITEM opm-core_SOURCES
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/LinearSolverIstl.cpp
		)
endif (NOT dune-istl_FOUND)
if (NOT SuiteSparse_FOUND)
	list (REMOVE_ITEM opm-core_SOURCES
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/call_umfpack.c
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/LinearSolverUmfpack.cpp
		)
	list (REMOVE_ITEM examples_SOURCES
		${PROJECT_SOURCE_DIR}/tutorials/tutorial2.cpp
		${PROJECT_SOURCE_DIR}/tutorials/tutorial3.cpp
		${PROJECT_SOURCE_DIR}/tutorials/tutorial4.cpp
		)
endif (NOT SuiteSparse_FOUND)

# we are not supposed to include the TinyXML test prog. regardless
list (REMOVE_ITEM opm-core_SOURCES
	${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/utility/parameters/tinyxml/xmltest.cpp
	)

# remove inline TinyXML if a system version was found
if (TinyXML_FOUND)
	file (GLOB_RECURSE _inline_tinyxml "${opm-core_DIR}/core/utility/parameters/tinyxml/*")
	foreach (_file IN LISTS _inline_tinyxml)
		list (REMOVE_ITEM opm-core_SOURCES ${_file})
	endforeach (_file)
endif (TinyXML_FOUND)

# anyhow remove it from the header list (so it doesn't get installed)
list (REMOVE_ITEM opm-core_HEADERS "${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/utility/parameters/tinyxml/tinystr.h")
list (REMOVE_ITEM opm-core_HEADERS "${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/utility/parameters/tinyxml/tinyxml.h")

# HAVE_ERT is used as an #ifdef, not as an #if in the source code, if it
# is not true, then it should be unset altogether
if (NOT HAVE_ERT)
	set (HAVE_ERT)
	list (REMOVE_ITEM examples_SOURCES
		${PROJECT_SOURCE_DIR}/examples/import_rewrite.cpp
		)
endif (NOT HAVE_ERT)
### --- end opm-core specific --- ###

# create configuration header which describes available features
# necessary to compile this library. singular version is the names that
# is required by this project alone, plural version transitively
# includes the necessary defines by the dependencies
include (ConfigVars)
list (APPEND ${project}_CONFIG_VARS ${${project}_CONFIG_VAR})

# write configuration variables to this file. note that it is a temporary.
message (STATUS "Writing config file \"${PROJECT_BINARY_DIR}/config.h\"...")
set (CONFIG_H "${PROJECT_BINARY_DIR}/config.h.tmp")
configure_vars (
	FILE  CXX  ${CONFIG_H}
	WRITE ${${project}_CONFIG_VARS}
	)

include (UseFortranWrappers)
define_fc_func (
	APPEND ${CONFIG_H}
	IF HAVE_AGMG # HAVE_BLAS HAVE_LAPACK
	)

# overwrite the config.h that is used by the code only if we have some
# real changes. thus, we don't have to recompile if a reconfigure is run
# due to some files being added, for instance
execute_process (COMMAND
	${CMAKE_COMMAND} -E copy_if_different ${CONFIG_H} ${PROJECT_BINARY_DIR}/config.h
	)

# compile main library; pull in all required includes and libraries
include (OpmCompile)
opm_compile (${project})

# installation target: copy the library together with debug and
# configuration files to system directories
include (OpmInstall)
opm_install (${project})
message (STATUS "This build defaults to installing in ${CMAKE_INSTALL_PREFIX}")

# installation of CMake modules to help user programs locate the library
include (OpmProject)
opm_cmake_config (${project})

# routines to build satellites such as tests, tutorials and samples
include (OpmSatellites)

# example programs are found in the tutorials/ and examples/ directory
opm_compile_satellites (${project} examples "" "")

# infrastructure for testing
enable_testing ()
include (CTest)

### --- begin opm-core specific --- ###
# conditionally disable tests when features aren't available
macro (cond_disable_test name)
	if ((NOT DEFINED HAVE_${name}) OR (NOT HAVE_${name}))
		message (STATUS "${name} test disabled, since ${name} is not found.")
		string (TOLOWER "${name}" name_lower)
		get_filename_component (test_${name}_FILE "tests/test_${name_lower}.cpp" ABSOLUTE)
		list (REMOVE_ITEM tests_SOURCES "${test_${name}_FILE}")
	endif ((NOT DEFINED HAVE_${name}) OR (NOT HAVE_${name}))
endmacro (cond_disable_test name)
cond_disable_test ("AGMG")
cond_disable_test ("ERT")
### --- end opm-core specific --- ###

# make datafiles necessary for tests available in output directory
opm_data (tests datafiles "${tests_DIR}")
opm_compile_satellites (${project} tests "" "${tests_REGEXP}")

# use this target to run all tests
add_custom_target (check
	COMMAND ${CMAKE_CTEST_COMMAND}
	DEPENDS tests
	COMMENT "Checking if library is functional"
	VERBATIM
	)

# generate documentation from source code with Doxygen;
# setup install target for this documentation
include (OpmDoc)
opm_doc (${project} ${doxy_dir})

# provide compatibility with using this build in dunecontrol
include (DuneCompat)
include (LibtoolArchives)
configure_la (${project} ${${project}_TARGET} ${project}_LIBTOOL_ARCHIVE)

### clean in-source builds ###
include (OpmDistClean)
opm_dist_clean (${project})
