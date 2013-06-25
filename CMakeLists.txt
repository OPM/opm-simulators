# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# Key information about the library
set (project "opm-autodiff")
set (${project}_NAME "${project}")
set (${project}_DESCRIPTION "OPM module for automatic differentiation")
set (${project}_DIR "opm")
set (${project}_VERSION_MAJOR 1)
set (${project}_VERSION_MINOR 0)
set (doxy_dir "doc/doxygen")

# Defines that must be present in config.h for our headers
set (${project}_CONFIG_VAR
	HAVE_DYNAMIC_BOOST_TEST
	)

# Prerequisites
set (${project}_DEPS
	# Compile with C99 support if available
	"C99"
	# Compile with C++0x/11 support if available
	"CXX11Features"
	# Various runtime library enhancements
	"Boost 1.39.0
		COMPONENTS system unit_test_framework REQUIRED"
	# DUNE prerequisites
	"dune-common REQUIRED;
	dune-istl REQUIRED;
	opm-core REQUIRED"
	# Eigen
	"Eigen3 3.1 REQUIRED"
	)

# Additional search modules
set (${project}_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
list (APPEND CMAKE_MODULE_PATH ${${project}_MODULE_DIR})

# Include special
if (CMAKE_VERSION VERSION_LESS "2.8.5")
	message (STATUS "Enabling compatibility modules for CMake 2.8.5")
	list (APPEND CMAKE_MODULE_PATH "${${project}_MODULE_DIR}/compat-2.8.5")
endif (CMAKE_VERSION VERSION_LESS "2.8.5")

if (CMAKE_VERSION VERSION_LESS "2.8.7")
	message (STATUS "Enabling compatibility modules for CMake 2.8.7")
	list (APPEND CMAKE_MODULE_PATH "${${project}_MODULE_DIR}/compat-2.8.7")
endif (CMAKE_VERSION VERSION_LESS "2.8.7")

# Don't write default flags into the cache, preserve that for user set values
include (AddOptions)
no_default_options ()

# C++ project
cmake_minimum_required (VERSION 2.8)
project (${${project}_NAME})
enable_language (C)
enable_language (CXX)

# Print system information to better pinpoint issues from log alone
include (UseSystemInfo)
system_info ()

# Very early try to print repo id (to pinpoint version if something goes wrong)
include (UseVCSInfo)
vcs_info ()

# Print toolchain information to identify compilers with potential bugs
include (UseCompVer)
compiler_info ()

# Default settings: build static debug library
include (OpmDefaults)
opm_defaults (${project})
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Use tricks to do faster builds
include (UseFastBuilds)

# Precompiled headers
include (UsePrecompHeaders)

# Optimise if we're not doing a debug build
include (UseOptimization)

# Turn on all warnings; this must be done before adding any
# prerequisites, in case they alter the list of warnings
include (UseWarnings)

# Parallel computing must be explicitly enabled
option (USE_MPI "Use Message Passing Interface for parallel computing" OFF)
if (NOT USE_MPI)
	set (CMAKE_DISABLE_FIND_PACKAGE_MPI TRUE)
endif (NOT USE_MPI)

### --- begin opm-autodiff specific --- ###
# Parallel programming
include (UseOpenMP)
find_openmp (${project})
### --- end opm-autodiff specific --- ###

# Macro to set standard variables (INCLUDE_DIRS, LIBRARIES etc.)
include (OpmFind)
find_and_append_package_list_to (${project} ${${project}_DEPS})

# Remove testing framework requirement from the main library;
# It is not possible to query for Boost twice with different components.
list (REMOVE_ITEM ${project}_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Don't import more libraries than we need to
include (UseOnlyNeeded)

# Put debug information into every executable
include (UseDebugSymbols)

# Detect if Boost is in a shared library
include (UseDynamicBoost)

# Needed for Debian installation scheme
include (UseMultiArch)

# This module contains code to figure out which files is where
include (OpmFiles)
opm_auto_dirs ()

# Put libraries in lib/
opm_out_dirs ()

# Identify the compilation units in the library; sources in opm/,
# tests files in tests/, examples in examples/
opm_sources (${project})

# Create configuration header which describes available features
# necessary to compile this library. Singular version is the names that
# is required by this project alone, plural version transitively
# includes the necessary defines by the dependencies
include (ConfigVars)
list (APPEND ${project}_CONFIG_VARS ${${project}_CONFIG_VAR})

# Write configuration variables to this file. Note that it is a temporary.
message (STATUS "Writing config file \"${PROJECT_BINARY_DIR}/config.h\"...")
set (CONFIG_H "${PROJECT_BINARY_DIR}/config.h.tmp")
configure_vars (
	FILE  CXX  ${CONFIG_H}
	WRITE ${${project}_CONFIG_VARS}
	)

# Overwrite the config.h that is used by the code only if we have some
# real changes. Thus, we don't have to recompile if a reconfigure is run
# due to some files being added, for instance
execute_process (COMMAND
	${CMAKE_COMMAND} -E copy_if_different ${CONFIG_H} ${PROJECT_BINARY_DIR}/config.h
	)

# Compile main library; pull in all required includes and libraries
include (OpmCompile)
opm_compile (${project})

# Installation target: copy the library together with debug and
# configuration files to system directories
include (OpmInstall)
opm_install (${project})
message (STATUS "This build defaults to installing in ${CMAKE_INSTALL_PREFIX}")

# Installation of CMake modules to help user programs locate the library
include (OpmProject)
opm_cmake_config (${project})

# Routines to build satellites such as tests, tutorials and samples
include (OpmSatellites)

# Example programs are found in the examples/ directory
opm_compile_satellites (${project} examples "" "")

# Infrastructure for testing
enable_testing ()
include (CTest)

# Make datafiles necessary for tests available in output directory
opm_data (tests datafiles "${tests_DIR}")
opm_compile_satellites (${project} tests "" "${tests_REGEXP}")

# Use this target to run all tests
add_custom_target (check
	COMMAND ${CMAKE_CTEST_COMMAND}
	DEPENDS tests
	COMMENT "Checking if library is functional"
	VERBATIM
	)

# Generate documentation from source code with Doxygen;
# setup install target for this documentation
include (OpmDoc)
opm_doc (${project} ${doxy_dir})

# Provide compatibility with using this build in dunecontrol
include (DuneCompat)
include (LibtoolArchives)
configure_la (${project} ${${project}_TARGET} ${project}_LIBTOOL_ARCHIVE)

### Clean in-source builds ###
include (OpmDistClean)
opm_dist_clean (${project})
