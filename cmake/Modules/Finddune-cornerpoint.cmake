# - Find OPM corner-point grid library
#
# Defines the following variables:
#   dune-cornerpoint_INCLUDE_DIRS    Directory of header files
#   dune-cornerpoint_LIBRARIES       Directory of shared object files
#   dune-cornerpoint_DEFINITIONS     Defines that must be set to compile
#   dune-cornerpoint_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_DUNE_CORNERPOINT            Binary value to use in config.h

# Copyright (C) 2013 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

include (dune-cornerpoint-prereqs)
include (OpmPackage)
find_opm_package (
  # module name
  "dune-cornerpoint"

  # dependencies
  "${dune-cornerpoint_DEPS}"
  
  # header to search for
  "dune/grid/CpGrid.hpp"

  # library to search for
  "dunecornerpoint"

  # defines to be added to compilations
  ""

  # test program
"#include <dune/grid/CpGrid.hpp>
int main (void) {
  Dune::CpGrid g;
  return 0;
}
"
  # config variables
  "${dune-cornerpoint_CONFIG_VAR}"
  )

#debug_find_vars ("dune-cornerpoint")
