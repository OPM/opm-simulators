# - Find OPM automatic differentiation library
#
# Defines the following variables:
#   opm-autodiff_INCLUDE_DIRS    Directory of header files
#   opm-autodiff_LIBRARIES       Directory of shared object files
#   opm-autodiff_DEFINITIONS     Defines that must be set to compile
#   opm-autodiff_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_OPM_AUTODIFF            Binary value to use in config.h

# Copyright (C) 2012 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

# use the generic find routine
include (opm-autodiff-prereqs)
include (OpmPackage)
find_opm_package (
  # module name
  "opm-autodiff"

  # dependencies
  "${opm-autodiff_DEPS}"
  
  # header to search for
  "opm/autodiff/AutoDiff.hpp"

  # library to search for
  "opmautodiff"

  # defines to be added to compilations
  ""

  # test program
"#include <opm/autodiff/AutoDiff.hpp>
int main (void) {
  Opm::AutoDiff<double> x = Opm::AutoDiff<double>::constant(42.);
  (void) x;
  return 0;  
}
"
  # config variables
  "${opm-autodiff_CONFIG_VAR}"
  )
include (UseDynamicBoost)
#debug_find_vars ("opm-autodiff")
