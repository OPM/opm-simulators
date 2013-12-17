# - Find OPM autodiff grid library
#
# Defines the following variables:
#   opm-autodiff_INCLUDE_DIRS    Directory of header files
#   opm-autodiff_LIBRARIES       Directory of shared object files
#   opm-autodiff_DEFINITIONS     Defines that must be set to compile
#   opm-autodiff_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_OPM_AUTODIFF            Binary value to use in config.h

# Copyright (C) 2013 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

include (opm-autodiff-prereqs)
include (OpmPackage)
find_opm_package (
  # module name
  "opm-autodiff"

  # dependencies
  "${opm-autodiff_DEPS}"

  # header to search for
  "opm/autodiff/SinglePhaseUpscaler.hpp"

  # library to search for
  "opmautodiff"

  # defines to be added to compilations
  ""

  # test program
"#include <opm/autodiff/AutoDiffBlock.hpp>
int main (void) {
  return 0;
}
"
  # config variables
  "${opm-autodiff_CONFIG_VAR}"
  )

#debug_find_vars ("opm-autodiff")
