# - Find routine for OPM-like modules
#
# Synopsis:
#
#   find_opm_package (module deps header lib defs prog conf)
#
# where
#
#   module    Name of the module, e.g. "dune-common"; this will be the
#             stem of all variables defined (see below).
#   deps      Semi-colon-separated list of dependent modules which must
#             be present; those that are required must be marked as such
#	          explicitly. Quote if more than one word is necessary to
#	          describe the dependency.
#   header    Name of the header file to probe for, e.g.
#             "dune/common/fvector.hh". Note that you should have to same
#             relative path here as is used in the header files.
#   lib       Name of the library to probe for, e.g. "dunecommon"
#   defs      Symbols that should be passed to compilations
#   prog      Program that should compile if library is present
#   conf      Symbols that should be present in config.h
#
# It will provide these standard Find-module variables:
#
#   ${module}_INCLUDE_DIRS    Directory of header files
#   ${module}_LIBRARIES       Directory of shared object files
#   ${module}_DEFINITIONS     Defines that must be set to compile
#   ${module}_CONFIG_VARS     List of defines that should be in config.h
#	${module}_QUIET           Verbosity of last find of this module
#   HAVE_${MODULE}            Binary value to use in config.h
#
# Note: Arguments should be quoted, otherwise a list will spill into the
#       next argument!

# Copyright (C) 2012 Uni Research AS
# This file is licensed under the GNU General Public License v3.0

# <http://www.vtk.org/Wiki/CMake:How_To_Find_Libraries>

include (Duplicates)
# append all items from src into dst; both must be *names* of lists
macro (append_found src dst)
  foreach (_item IN LISTS ${src})
	if (NOT "${_item}" MATCHES "-NOTFOUND$")
	  list (APPEND ${dst} ${_item})
	endif (NOT "${_item}" MATCHES "-NOTFOUND$")
  endforeach (_item)
endmacro (append_found src dst)

function (find_opm_package module deps header lib defs prog conf)
  # variables to pass on to other packages
  if (FIND_QUIETLY)
	set (${module}_QUIET "QUIET")
  else (FIND_QUIETLY)
	set (${module}_QUIET "")
  endif (FIND_QUIETLY)

  # if someone else has included this test, don't do it again
  if (${${module}_FOUND})
	return ()
  endif (${${module}_FOUND})

  # dependencies on other packages; underscore version only contains the
  # name of the other package
  set (_deps)
  foreach (_dep IN LISTS deps)
	separate_arguments (_args UNIX_COMMAND ${_dep})
	find_package (${_args} QUIET)
	list (GET _args 0 _name_only)
	list (APPEND _deps ${_name_only})
  endforeach (_dep)

  # see if there is a pkg-config entry for this package, and use those
  # settings as a starting point
  find_package (PkgConfig)
  pkg_check_modules (PkgConf_${module} QUIET ${module})

  # these variables have non-standard names in FindPkgConfig (sic)
  set (${module}_DEFINITIONS ${PkgConf_${module}_CFLAGS_OTHER})
  set (${module}_LINKER_FLAG ${PkgConf_${module}_LDFLAGS_OTHER})

  # in addition to accepting mod-ule_ROOT, we also accept the somewhat
  # more idiomatic MOD_ULE_ROOT variant
  string (TOUPPER ${module} MODULE_UPPER)
  string (REPLACE "-" "_" MODULE ${MODULE_UPPER})

  # if the user hasn't specified any location, and it isn't found
  # in standard system locations either, then start to wander
  # about and look for it in proximity to ourself. Qt Creator likes
  # to put the build-directories as siblings to the source trees,
  # but with a -build suffix
  if (NOT (${module}_DIR OR ${module}_ROOT OR ${MODULE}_ROOT))
	string (TOLOWER "${module}" _module_lower)
	set (_guess
	  "../${module}"
	  "../${module}-build"
	  "../${_module_lower}"
	  "../${_module_lower}-build"
	  )
	set (_guess_bin)
	foreach (_item IN ITEMS ${_guess})
	  list (APPEND _guess_bin "${PROJECT_BINARY_DIR}/${_item}")
	endforeach (_item)
  endif (NOT (${module}_DIR OR ${module}_ROOT OR ${MODULE}_ROOT))

  # search for this include and library file to get the installation
  # directory of the package; hints are searched before the system locations,
  # paths are searched afterwards
  find_path (${module}_INCLUDE_DIR
	NAMES "${header}"
	PATHS ${_guess}
	HINTS ${${module}_DIR} ${${module}_ROOT} ${${MODULE}_ROOT} ${PkgConf_${module}_INCLUDE_DIRS}
	PATH_SUFFIXES "include"
	)

  # some modules are all in headers
  if (NOT "${lib}" STREQUAL "")
	find_library (${module}_LIBRARY
	  NAMES "${lib}"
	  PATHS ${_guess_bin}
	  HINTS ${${module}_DIR} ${${module}_ROOT} ${${MODULE}_ROOT} ${PkgConf_${module}_LIBRARY_DIRS}
	  PATH_SUFFIXES "lib" "lib/.libs" ".libs" "lib32" "lib64" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
	  )
  else (NOT "${lib}" STREQUAL "")
	set (${module}_LIBRARY "")
  endif (NOT "${lib}" STREQUAL "")

  # add dependencies so that our result variables are complete
  # list of necessities to build with the software
  set (${module}_INCLUDE_DIRS "${${module}_INCLUDE_DIR}")
  set (${module}_LIBRARIES "${${module}_LIBRARY}")
  foreach (_dep IN LISTS _deps)
	# only add those packages we actually found (find_package will show
	# an error if it was marked as REQUIRED)
	if (${_dep}_FOUND)
	  list (APPEND ${module}_INCLUDE_DIRS ${${_dep}_INCLUDE_DIRS})
	  list (APPEND ${module}_LIBRARIES ${${_dep}_LIBRARIES})
	  list (APPEND ${module}_DEFINITIONS ${${_dep}_DEFINITIONS})
	  list (APPEND ${module}_CONFIG_VARS ${${_dep}_CONFIG_VARS})
	  list (APPEND ${module}_LINKER_FLAGS ${${_dep}_LINKER_FLAGS})
	endif (${_dep}_FOUND)
  endforeach (_dep)

  # compile with this option to avoid avalanche of warnings
  set (${module}_DEFINITIONS "${${module}_DEFINITIONS}")
  foreach (_def IN LISTS defs)
	list (APPEND ${module}_DEFINITIONS "-D${_def}")
  endforeach (_def)

  # tidy the lists before returning them
  list (REMOVE_DUPLICATES ${module}_INCLUDE_DIRS)
  remove_duplicate_libraries (${module})
  list (REMOVE_DUPLICATES ${module}_DEFINITIONS)

  # these defines are used in dune/${module} headers, and should be put
  # in config.h when we include those
  foreach (_var IN LISTS conf)
	# massage the name to remove source code formatting
	string (REGEX REPLACE "^[\n\t\ ]+" "" _var "${_var}")
	string (REGEX REPLACE "[\n\t\ ]+$" "" _var "${_var}")
	list (APPEND ${module}_CONFIG_VARS ${_var})
  endforeach (_var)
  foreach (_dep in _deps)
	if (DEFINED ${_dep}_CONFIG_VARS)
	  list (APPEND ${module}_CONFIG_VARS ${_dep}_CONFIG_VARS)
	endif (DEFINED ${_dep}_CONFIG_VARS)
  endforeach (_dep)
  if (${module}_CONFIG_VARS)
	list (REMOVE_DUPLICATES ${module}_CONFIG_VARS)
  endif (${module}_CONFIG_VARS)

  # these are the defines that should be set when compiling
  # without config.h
  config_cmd_line (${module}_CMD_CONFIG ${module}_CONFIG_VARS)

  # check that we can compile a small test-program
  include (CMakePushCheckState)
  cmake_push_check_state ()
  include (CheckCXXSourceCompiles)
  # only add these if they are actually found; otherwise it won't
  # compile and the variable won't be set
  append_found (${module}_INCLUDE_DIRS CMAKE_REQUIRED_INCLUDES)
  append_found (${module}_LIBRARIES CMAKE_REQUIRED_LIBRARIES)
  # since we don't have any config.h yet
  list (APPEND CMAKE_REQUIRED_DEFINITIONS ${${module}_DEFINITIONS})
  list (APPEND CMAKE_REQUIRED_DEFINITIONS ${${module}_CMD_CONFIG})
  check_cxx_source_compiles ("${prog}" HAVE_${MODULE})
  cmake_pop_check_state ()

  # write status message in the same manner as everyone else
  include (FindPackageHandleStandardArgs)
  if ("${lib}" STREQUAL "")
	set (_lib_var "")
  else ("${lib}" STREQUAL "")
	set (_lib_var "${module}_LIBRARY")
  endif ("${lib}" STREQUAL "")
  find_package_handle_standard_args (
	${module}
	DEFAULT_MSG
	${module}_INCLUDE_DIR ${_lib_var} HAVE_${MODULE}
	)

  # allow the user to override these from user interface
  mark_as_advanced (${module}_INCLUDE_DIR)
  mark_as_advanced (${module}_LIBRARY)

  # some genius that coded the FindPackageHandleStandardArgs figured out
  # that the module name should be in uppercase (?!)
  set (${module}_FOUND "${${MODULE_UPPER}_FOUND}" PARENT_SCOPE)

  # return these variables to the caller
  set (${module}_INCLUDE_DIRS "${${module}_INCLUDE_DIRS}" PARENT_SCOPE)
  set (${module}_LIBRARIES "${${module}_LIBRARIES}" PARENT_SCOPE)
  set (${module}_DEFINITIONS "${${module}_DEFINITIONS}" PARENT_SCOPE)
  set (${module}_CONFIG_VARS "${${module}_CONFIG_VARS}" PARENT_SCOPE)
  set (${module}_LINKER_FLAGS "${${module}_LINKER_FLAGS}" PARENT_SCOPE)
  set (${module}_QUIET "${${module}_QUIET}" PARENT_SCOPE)
  set (HAVE_${MODULE} "${HAVE_${MODULE}}" PARENT_SCOPE)
endfunction (find_opm_package module deps header lib defs prog conf)

# print all variables defined by the above macro
function (debug_find_vars module)
  message (STATUS "${module}_FOUND        = ${${module}_FOUND}")
  message (STATUS "${module}_INCLUDE_DIRS = ${${module}_INCLUDE_DIRS}")
  message (STATUS "${module}_LIBRARIES    = ${${module}_LIBRARIES}")
  message (STATUS "${module}_DEFINITIONS  = ${${module}_DEFINITIONS}")
  message (STATUS "${module}_CONFIG_VARS  = ${${module}_CONFIG_VARS}")
  message (STATUS "${module}_LINKER_FLAGS = ${${module}_LINKER_FLAGS}")
  message (STATUS "${module}_QUIET        = ${${module}_QUIET}")
  string (TOUPPER ${module} MODULE)
  string (REPLACE "-" "_" MODULE ${MODULE})  
  message (STATUS "HAVE_${MODULE}         = ${HAVE_${MODULE}}")
endfunction (debug_find_vars module)

# generate a command-line that can be used to pass variables before
# config.h is available (such as probe tests). varname is the *name*
# of the variable to receive the result, defs is a list of the *names*
# which should be passed
function (config_cmd_line varname defs)
  # process each variable
  foreach (_var IN LISTS ${defs})
	# only generate an entry if the define was actually set
	if ((DEFINED ${_var}) AND (NOT "${${_var}}" STREQUAL ""))
	  # numbers are not quoted, strings are
	  if (${_var} MATCHES "[0-9]+")
		set (_quoted "${${_var}}")
	  else (${_var} MATCHES "[0-9]+")
		set (_quoted "\"${${_var}}\"")
	  endif (${_var} MATCHES "[0-9]+")
	  # add command-line option to define this variable
	  list (APPEND _cmdline "-D${_var}=${_quoted}")
	endif ((DEFINED ${_var}) AND (NOT "${${_var}}" STREQUAL ""))
  endforeach (_var)
  # return the resulting command-line options for defining vars
  set (${varname} "${_cmdline}" PARENT_SCOPE)
endfunction (config_cmd_line)
