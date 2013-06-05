# - Compile main library target

macro (opm_compile opm)
  # some CMake properties do not do list expansion
  string (REPLACE ";" " " ${opm}_LINKER_FLAGS_STR "${${opm}_LINKER_FLAGS}")

  # name of the library should not contain dashes, as CMake will
  # define a symbol with that name, and those cannot contain dashes
  string (REPLACE "-" "" ${opm}_TARGET "${${opm}_NAME}")

  # all public header files are together with the source. prepend our own
  # source path to the one of the dependencies so that our version of any
  # ambigious paths are used.
  set (${opm}_INCLUDE_DIR "${PROJECT_SOURCE_DIR}")
  set (${opm}_INCLUDE_DIRS ${${opm}_INCLUDE_DIR} ${${opm}_INCLUDE_DIRS})

  # create this library
  include_directories (${${opm}_INCLUDE_DIRS})
  link_directories (${${opm}_LIBRARY_DIRS})
  add_definitions (${${opm}_DEFINITIONS})
  add_library (${${opm}_TARGET} ${${opm}_LIBRARY_TYPE} ${${opm}_SOURCES})
  set (${opm}_VERSION "${${opm}_VERSION_MAJOR}.${${opm}_VERSION_MINOR}")
  set_target_properties (${${opm}_TARGET} PROPERTIES
	SOVERSION ${${opm}_VERSION_MAJOR}
	VERSION ${${opm}_VERSION}
	LINK_FLAGS "${${opm}_LINKER_FLAGS_STR}"
	)
  target_link_libraries (${${opm}_TARGET} ${${opm}_LIBRARIES})

  # queue this executable to be stripped
  strip_debug_symbols (${${opm}_TARGET} ${opm}_DEBUG)

  # pre-compile common headers; this is setup *after* the library to pick
  # up extra options set there
  if (PRECOMPILE_HEADERS)
	get_target_property (_type ${${opm}_TARGET} TYPE)
	precompile_header (CXX ${_type}
	  HEADER "${${opm}_PRECOMP_CXX_HEADER}"
	  TARGET ${opm}_CXX_pch
	  FLAGS  ${opm}_PRECOMP_CXX_FLAGS
	  )
	# must set property on source files instead of entire target, because
	# it only applies to C++ modules (and cannot be used for C)
	set_source_files_properties (${${opm}_CXX_SOURCES} PROPERTIES
	  OBJECT_DEPENDS "${${opm}_CXX_pch}"
	  COMPILE_FLAGS  "${${opm}_PRECOMP_CXX_FLAGS}"
	  )
	message (STATUS "Precompiled headers: ${${opm}_CXX_pch}")
  endif (PRECOMPILE_HEADERS)

  # we need to know the name of the library which is generated
  get_target_property (${opm}_LIBRARY ${${opm}_TARGET} LOCATION)
  
endmacro (opm_compile opm)
