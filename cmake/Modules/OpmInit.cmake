# - Initialize project-specific variables
#
# This will read the dune.module file for project information and
# set the following variables:
#
#	project                     From the Module: field
#	${project}_NAME             Same as above
#	${project}_DESCRIPTION      From the Description: field
#	${project}_VERSION_MAJOR    From the Version: field
#	${project}_VERSION_MINOR    From the Version: field also
#
# This module should be the first to be included in the project,
# because most of the others (OpmXxx.cmake) use these variables.

# helper macro to retrieve a single field of a dune.module file
macro(OpmGetDuneModuleDirective field variable contents)
  string (REGEX MATCH ".*${field}:[ ]*([^\n]+).*" ${variable} "${contents}")
  string (REGEX REPLACE ".*${field}:[ ]*([^\n]+).*" "\\1" "${variable}" "${${variable}}")
  string (STRIP "${${variable}}" ${variable})
endmacro()

function (OpmInitProjVars)
  # locate the "dune.module" file
  set (DUNE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/dune.module")

  # read this file into a variable
  file (READ "${DUNE_MODULE_PATH}" DUNE_MODULE)

  # read fields from the file
  OpmGetDuneModuleDirective ("Module" project "${DUNE_MODULE}")
  OpmGetDuneModuleDirective ("Description" description "${DUNE_MODULE}")
  OpmGetDuneModuleDirective ("Version" version "${DUNE_MODULE}")

  # parse the version number
  set (verno_regex "^([0-9]*)\\.([0-9]*).*\$")
  string (REGEX REPLACE "${verno_regex}" "\\1" major "${version}")
  string (REGEX REPLACE "${verno_regex}" "\\2" minor "${version}")

  # return these variables
  set (project "${project}" PARENT_SCOPE)
  set (${project}_NAME "${project}" PARENT_SCOPE)
  set (${project}_DESCRIPTION "${description}" PARENT_SCOPE)
  set (${project}_VERSION_MAJOR "${major}" PARENT_SCOPE)
  set (${project}_VERSION_MINOR "${minor}" PARENT_SCOPE)  
endfunction ()

macro (OpmInitDirVars)
  # these are the most common (and desired locations)
  set (${project}_DIR "opm")
  set (doxy_dir "doc/doxygen")

  # but for backward compatibility we can override it
  if (COMMAND dir_hook)
	dir_hook ()	
  endif (COMMAND dir_hook)
endmacro ()

OpmInitProjVars ()
OpmInitDirVars ()
