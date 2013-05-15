# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# key information about the library
set (project "opm-autodiff")
set (${project}_NAME "${project}")
set (${project}_DESCRIPTION "OPM module for automatic differentiation")
set (${project}_DIR "opm")
set (${project}_VERSION_MAJOR 1)
set (${project}_VERSION_MINOR 0)
set (doxy_dir "doc/doxygen")

# defines that must be present in config.h for our headers
set (${project}_CONFIG_VAR
	HAVE_DYNAMIC_BOOST_TEST
	)

# dependencies
set (${project}_DEPS
	# compile with C99 support if available
	"C99"
	# compile with C++0x/11 support if available
	"CXX11Features"
	# various runtime library enhancements
	"Boost 1.39.0
		COMPONENTS system unit_test_framework REQUIRED"
	# DUNE dependency
	"dune-common REQUIRED;
	dune-istl REQUIRED;
	opm-core REQUIRED"
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

### --- begin opm-porsol specific --- ###
# parallel programming
include (UseOpenMP)
find_openmp (${project})
### --- end opm-porsol specific --- ###

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
# tests files in tests/, examples in examples/
opm_sources (${project})

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

# example programs are found in the examples/ directory
opm_compile_satellites (${project} examples "" "")

# infrastructure for testing
enable_testing ()
include (CTest)

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
