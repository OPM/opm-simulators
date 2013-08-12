# - Find DUNE localfunctions library
#
# Defines the following variables:
#   dune-localfunctions_INCLUDE_DIRS      Directory of header files
#   dune-localfunctions_LIBRARIES         Directory of shared object files
#   dune-localfunctions_DEFINITIONS       Defines that must be set to compile
#   dune-localfunctions_CONFIG_VARS       List of defines that should be in config.h
#   HAVE_DUNE_LOCALFUNCTIONS              Binary value to use in config.h

# Copyright (C) 2012 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

# dune-common is only required if dune-localfunctions is; the "required-ness" is
# not transitive as far as CMake is concerned (i.e. if an optional package
# requests a package to be required, the build will fail if it's not found)
if (dune-localfunctions_FIND_REQUIRED)
  set (_require_dune_common "REQUIRED")
endif (dune-localfunctions_FIND_REQUIRED)

include (OpmPackage)
find_opm_package (
  # module name
  "dune-localfunctions"

  # required dependencies
  "dune-common ${_require_dune_common}"
  # header to search for
  "dune/localfunctions/common/localbasis.hh"

  # library to search for
  ""

  # defines to be added to compilations
  ""

  # test program
"#include   <dune/localfunctions/common/localbasis.hh>

int main (void) {
  return 0;
}
"
  # config variables
  "")
#debug_find_vars ("dune-localfunctions")

# make version number available in config.h
include (UseDuneVer)
find_dune_version ("dune" "localfunctions")
