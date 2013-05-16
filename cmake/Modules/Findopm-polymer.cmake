# - Find OPM polymer library
#
# Defines the following variables:
#   opm-polymer_INCLUDE_DIRS    Directory of header files
#   opm-polymer_LIBRARIES       Directory of shared object files
#   opm-polymer_DEFINITIONS     Defines that must be set to compile
#   opm-polymer_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_OPM_POLYMER            Binary value to use in config.h

# Copyright (C) 2013 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

include (OpmPackage)
find_opm_package (
  # module name
  "opm-polymer"

  # dependencies
  "ERT;
  opm-core REQUIRED
  "
  # header to search for
  "opm/polymer/PolymerState.hpp"

  # library to search for
  "opmpolymer"

  # defines to be added to compilations
  ""

  # test program
"#include <opm/polymer/PolymerState.hpp>
int main (void) {
  Opm::PolymerState s;
  return 0;
}
"
  # config variables
  "HAVE_ERT")

#debug_find_vars ("opm-polymer")
