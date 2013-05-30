# - Use only needed imports from libraries
#
# Add the -Wl,--as-needed flag to the default linker flags on Linux
# in order to get only the minimal set of dependencies.

function (prepend var_name value)
  if (${var_name})
	set (${var_name} "${value} ${${var_name}}" PARENT_SCOPE)
  else (${var_name})
	set (${var_name} "${value}")
  endif (${var_name})
endfunction (prepend var_name value)

# only ELF shared objects can be underlinked, and only GNU will accept
# these parameters; otherwise just leave it to the defaults
if ((CMAKE_CXX_PLATFORM_ID STREQUAL "Linux") AND CMAKE_COMPILER_IS_GNUCC)
  # these are the modules whose probes will turn up incompatible
  # flags on some systems
  set (_maybe_underlinked
	SuiteSparse
	)
  # check if any modules actually reported problems (by setting an
  # appropriate linker flag)
  set (_underlinked FALSE)
  foreach (_module IN LISTS _maybe_underlinked)
	if ("${${_module}_LINKER_FLAGS}" MATCHES "-Wl,--no-as-needed")
	  set (_underlinked TRUE)
	endif ("${${_module}_LINKER_FLAGS}" MATCHES "-Wl,--no-as-needed")
  endforeach (_module)
  # if we didn't have any problems, then go ahead and add
  if (NOT _underlinked)
	prepend (CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
	prepend (CMAKE_MODULE_LINKER_FLAGS "-Wl,--as-needed")
	prepend (CMAKE_SHARED_LINKER_FLAGS "-Wl,--as-needed")
  endif (NOT _underlinked)
endif ((CMAKE_CXX_PLATFORM_ID STREQUAL "Linux")  AND CMAKE_COMPILER_IS_GNUCC)
