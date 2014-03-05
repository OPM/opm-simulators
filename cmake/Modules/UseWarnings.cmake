# - Turn on warnings when compiling

include (AddOptions)
include (UseCompVer)
is_compiler_gcc_compatible ()

if (CXX_COMPAT_GCC)
  # default warnings flags, if not set by user
  set_default_option (CXX _warn_flag "-Wall" "(^|\ )-W")
  if (_warn_flag)
	message (STATUS "All warnings enabled: ${_warn_flag}")
	add_options (ALL_LANGUAGES ALL_BUILDS "${_warn_flag}")
  endif (_warn_flag)
endif ()

option(SILENCE_DUNE_WARNINGS "Disable warnings from DUNE?" OFF)
if(SILENCE_DUNE_WARNINGS AND CXX_COMPAT_GCC)
  file(WRITE ${CMAKE_BINARY_DIR}/dune_disable_pragmas.h "
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored \"-Wshadow\"
#pragma GCC diagnostic ignored \"-Wunused-parameter\"
#pragma GCC diagnostic ignored \"-Wignored-qualifiers\"
#pragma GCC diagnostic ignored \"-Wmismatched-tags\"")
  file(WRITE ${CMAKE_BINARY_DIR}/dune_reenable_pragmas.h "#pragma GCC diagnostic pop")
else()
    file(WRITE ${CMAKE_BINARY_DIR}/dune_disable_pragmas.h "")
    file(WRITE ${CMAKE_BINARY_DIR}/dune_reenable_pragmas.h "")
endif()
