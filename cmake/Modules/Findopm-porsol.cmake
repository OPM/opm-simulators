# - Find OPM porous media solver library
#
# Defines the following variables:
#   opm-porsol_INCLUDE_DIRS    Directory of header files
#   opm-porsol_LIBRARIES       Directory of shared object files
#   opm-porsol_DEFINITIONS     Defines that must be set to compile
#   opm-porsol_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_OPM_PORSOL            Binary value to use in config.h

# Copyright (C) 2013 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

include (OpmPackage)
find_opm_package (
  # module name
  "opm-porsol"

  # dependencies
  "dune-common REQUIRED;
  dune-grid REQUIRED;
  dune-istl REQUIRED;
  opm-core REQUIRED;
  dune-cornerpoint REQUIRED
  "
  # header to search for
  "opm/porsol/mimetic/IncompFlowSolverHybrid.hpp"

  # library to search for
  "opmporsol"

  # defines to be added to compilations
  ""

  # test program
"#include <opm/porsol/mimetic/IncompFlowSolverHybrid.hpp>
int main (void) {
  return 0;
}
"
  # config variables
  "HAVE_VALGRIND")

#debug_find_vars ("opm-porsol")
