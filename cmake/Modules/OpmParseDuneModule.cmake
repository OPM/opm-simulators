# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:
include(CMakeParseArguments)

# helper macro to retrieve a single field of a dune.module file
macro(OpmGetDuneModuleDirective_ FieldName OutputVariable DuneModuleContents)
  string(REGEX MATCH ".*${FieldName}:[ ]*([^\n]+).*" ${OutputVariable} "${DuneModuleContents}")
  string(REGEX REPLACE ".*${FieldName}:[ ]*([^\n]+).*" "\\1" "${OutputVariable}" "${${OutputVariable}}")
endmacro()

# This macro parses the dune.module file of a Dune module.
#
# Usage:
#
# OpmParseDuneModuleInfo(DuneModuleName
#              [FILE_NAME "PathTo/dune.module"])
#
# This function sets the following variables:
#
# ${MODULE_NAME}_NAME: The name of the module which ought to be used when referencing it in texts
# ${MODULE_NAME}_DESCRIPTION: A textual description of the module
# ${MODULE_NAME}_DEPENDS: The dependencies of the module
# ${MODULE_NAME}_URL: The URL of the module's website
# ${MODULE_NAME}_MAINTAINER_NAME: The real name of the module maintainer(s)
# ${MODULE_NAME}_MAINTAINER_EMAIL: The e-mail address of the module maintainer(s)
# ${MODULE_NAME}_VERSION: version string of the dune module
# ${MODULE_NAME}_VERSION_MAJOR: major version of the dune module
# ${MODULE_NAME}_VERSION_MINOR: minor version of the dune module
# ${MODULE_NAME}_VERSION_REVISION: revision of the dune module
# ${MODULE_NAME}_CODENAME: code name for the version of the module
function(OpmParseDuneModule ModuleName)
  CMAKE_PARSE_ARGUMENTS(
    CURMOD # prefix
    "" # flags
    "FILE_NAME" # one value args
    "" # multi-value args
    ${ARGN})

  # create an uppercase, underscore version of the given module name
  string(TOUPPER "${ModuleName}" MODULE_NAME)
  string(REPLACE "-" "_" MODULE_NAME "${MODULE_NAME}")
  
  if (NOT CURMOD_FILE_NAME)
    set(CURMOD_FILE_NAME "${${MODULE_NAME}_DIR}/dune.module")
  endif()

  # read the dune.module file
  file(READ "${CURMOD_FILE_NAME}" DUNE_MODULE)

  # set the module name
  OpmGetDuneModuleDirective_("Module" ${MODULE_NAME}_NAME "${DUNE_MODULE}")

  # set the module description
  OpmGetDuneModuleDirective_("Description" ${MODULE_NAME}_DESCRIPTION "${DUNE_MODULE}")

  # set the dependencies
  OpmGetDuneModuleDirective_("Depends" ${MODULE_NAME}_DEPENDS "${DUNE_MODULE}")

  # set the URL of the module's website
  OpmGetDuneModuleDirective_("Url" ${MODULE_NAME}_URL "${DUNE_MODULE}")

  # set the project version. also split the version string into MAJOR.MINOR.REVISON
  OpmGetDuneModuleDirective_("Version" ${MODULE_NAME}_VERSION "${DUNE_MODULE}")

  string(REGEX REPLACE "^([0-9]*)\\..*\$" "\\1" ${MODULE_NAME}_VERSION_MAJOR "${${MODULE_NAME}_VERSION}")
  string(REGEX REPLACE "^[0-9]*\\.([0-9]*).*\$" "\\1" ${MODULE_NAME}_VERSION_MINOR "${${MODULE_NAME}_VERSION}")
  string(REGEX REPLACE "^[0-9]*\\.[0-9]*\\.([0-9]*).*\$" "\\1" ${MODULE_NAME}_VERSION_REVISION "${${MODULE_NAME}_VERSION}")

  # if the regular expression for the revision did not match, we use "0"
  # as the revision number. (we silently assume, that the regexps for
  # the major and minor version match.)
  if ("${${MODULE_NAME}_VERSION_REVISION}" STREQUAL "${${MODULE_NAME}_VERSION}")
    set(${MODULE_NAME}_VERSION_REVISION "0")
  endif()

  # set the maintainer email (the default Dune autotools build system
  # assumes that dune.module's 'Maintainer' field only contains the
  # email address of the maintainer. Using the format 'Maintainer:
  # Maintainer Name <maintainer@address.org>' makes the DUNE autotools
  # build system choke, so we introduce a new field 'MaintainerName'
  # which is ignored by the DUNE autotools build system.)
  OpmGetDuneModuleDirective_("MaintainerName" ${MODULE_NAME}_MAINTAINER_NAME "${DUNE_MODULE}")
  OpmGetDuneModuleDirective_("Maintainer" ${MODULE_NAME}_MAINTAINER_EMAIL "${DUNE_MODULE}")

  # find codename string
  OpmGetDuneModuleDirective_("Codename" ${MODULE_NAME}_CODENAME "${DUNE_MODULE}")

  ################
  # export all output variables
  ################

  # needed for some other OPM CMake modules: Variable without dashes
  # replaced by underscores and non-uppercase
  set(${ModuleName}_NAME "${${MODULE_NAME}_NAME}" PARENT_SCOPE)
  set(${MODULE_NAME}_NAME "${${MODULE_NAME}_NAME}" PARENT_SCOPE)
  set(${MODULE_NAME}_DESCRIPTION "${${MODULE_NAME}_DESCRIPTION}" PARENT_SCOPE)
  set(${MODULE_NAME}_DEPENDS "${${MODULE_NAME}_DEPENDS}" PARENT_SCOPE)
  set(${MODULE_NAME}_URL "${${MODULE_NAME}_URL}" PARENT_SCOPE)
  set(${MODULE_NAME}_VERSION "${${MODULE_NAME}_VERSION}" PARENT_SCOPE)
  set(${MODULE_NAME}_VERSION_MAJOR "${${MODULE_NAME}_VERSION_MAJOR}" PARENT_SCOPE)
  set(${MODULE_NAME}_VERSION_MINOR "${${MODULE_NAME}_VERSION_MINOR}" PARENT_SCOPE)
  set(${MODULE_NAME}_VERSION_REVISION "${${MODULE_NAME}_VERSION_REVISION}" PARENT_SCOPE)
  set(${MODULE_NAME}_MAINTAINER_NAME "${${MODULE_NAME}_MAINTAINER_NAME}" PARENT_SCOPE)
  set(${MODULE_NAME}_MAINTAINER_EMAIL "${${MODULE_NAME}_MAINTAINER_EMAIL}" PARENT_SCOPE)
  set(${MODULE_NAME}_CODENAME "${${MODULE_NAME}_CODENAME}" PARENT_SCOPE)
endfunction()
