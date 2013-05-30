# - Turn on optimizations based on build type

include(TestCXXAcceptsFlag)
include (AddOptions)

# if we are building a debug target, then disable all optimizations
# otherwise, turn them on. indicate to the code what we have done
# so it can turn on assertions etc.
if (CMAKE_COMPILER_IS_GNUCXX)
  # use these options for debug builds - no optimizations
  add_options (
	ALL_LANGUAGES
	"Debug"
	"-O0" "-DDEBUG"
	)

  # extra flags passed for optimization
  set (_opt_flags "")

  # link-time (a.k.a. global) optimizations
  check_cxx_accepts_flag ("-flto" HAVE_LINK_OPTS)
  if (HAVE_LINK_OPTS)
	list (APPEND _opt_flags "-flto")
  endif (HAVE_LINK_OPTS)

  # native instruction set tuning
  option (WITH_NATIVE "Use native instruction set" ON)
  if (WITH_NATIVE)
	check_cxx_accepts_flag ("-mtune=native" HAVE_MTUNE)
	if (HAVE_MTUNE)
	  list (APPEND _opt_flags "-mtune=native")
	endif (HAVE_MTUNE)
  endif (WITH_NATIVE)

  # use these options for release builds - full optimization
  add_options (
	ALL_LANGUAGES
	"Release;RelWithDebInfo;MinSizeRel"
	"-O3" "-DNDEBUG" ${_opt_flags}
	)
endif (CMAKE_COMPILER_IS_GNUCXX)
