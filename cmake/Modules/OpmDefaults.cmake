# - Default settings for the build

include (UseCompVer)
is_compiler_gcc_compatible ()
include(TestCXXAcceptsFlag)

macro (opm_defaults opm)
  # if we are installing a development version (default when checking out of
  # VCS), then remember which directories were used when configuring. package
  # distribution should disable this option.
  option (USE_RUNPATH "Embed original dependency paths in installed library" ON)
  if (USE_RUNPATH)
	if (CXX_COMPAT_GCC)
	  check_cxx_accepts_flag ("-Wl,--enable-new-dtags" HAVE_RUNPATH)
	  if (HAVE_RUNPATH)
		list (APPEND ${opm}_LINKER_FLAGS "-Wl,--enable-new-dtags")
	  endif (HAVE_RUNPATH)
	endif ()
	# set this to avoid CMake stripping it off again
	set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif (USE_RUNPATH)

  # build release by default
  if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Release")
  endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

  # default to building a static library, but let user override
  if (DEFINED BUILD_SHARED_LIBS)
	if (BUILD_SHARED_LIBS)
	  set (${opm}_LIBRARY_TYPE SHARED)
	else (BUILD_SHARED_LIBS)
	  set (${opm}_LIBRARY_TYPE STATIC)
	endif (BUILD_SHARED_LIBS)
  else (DEFINED BUILD_SHARED_LIBS)
	set (${opm}_LIBRARY_TYPE STATIC)
  endif (DEFINED BUILD_SHARED_LIBS)

  # precompile standard headers to speed up compilation
  # unfortunately, this functionality is buggy and tends to segfault at
  # least up to version 4.7.2, so it should be disabled by default there
  set (_precomp_def OFF)
  option (PRECOMPILE_HEADERS "Precompile common headers for speed." ${_precomp_def})
  mark_as_advanced (PRECOMPILE_HEADERS)
  if (NOT PRECOMPILE_HEADERS)
	message (STATUS "Precompiled headers: disabled")
  endif(NOT PRECOMPILE_HEADERS)

  # Use of OpenMP is considered experimental
  set (USE_OPENMP_DEFAULT OFF)

  # if we are on a system where CMake 2.6 is the default (Hi RHEL 6!),
  # the configuration files for Boost will trip up the library paths
  # (look for /usr/lib64/lib64/ in the log) when used with FindBoost
  # module bundled with CMake 2.8. this can be circumvented by turning
  # off config mode probing if we have not explicitly specified a
  # directory to look for it. for more details, see
  # <http://stackoverflow.com/questions/9948375/cmake-find-package-succeeds-but-returns-wrong-path>
  if (NOT BOOST_ROOT)
	set (Boost_NO_BOOST_CMAKE ON)
  endif (NOT BOOST_ROOT)
endmacro (opm_defaults opm)
