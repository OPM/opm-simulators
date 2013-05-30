# - Find OPM core library
#
# Defines the following variables:
#   opm-core_INCLUDE_DIRS    Directory of header files
#   opm-core_LIBRARIES       Directory of shared object files
#   opm-core_DEFINITIONS     Defines that must be set to compile
#   opm-core_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_OPM_CORE            Binary value to use in config.h

# Copyright (C) 2012 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

# use the generic find routine
include (OpmPackage)
find_opm_package (
  # module name
  "opm-core"

  # dependencies
  "C99;
  CXX11Features;
  Boost 1.39.0
    COMPONENTS date_time filesystem system unit_test_framework REQUIRED;
  BLAS REQUIRED;
  LAPACK REQUIRED;
  SuiteSparse COMPONENTS umfpack;
  SUPERLU;
  TinyXML;
  ERT;
  dune-istl
  "
  # header to search for
  "opm/core/grid.h"

  # library to search for
  "opmcore"

  # defines to be added to compilations
  ""

  # test program
"#include <opm/core/grid.h>
int main (void) {
  struct UnstructuredGrid *g;
  g = create_grid_empty ();
  destroy_grid (g);
  return 0;  
}
"
  # config variables
  "HAVE_AGMG;
  HAVE_DUNE_ISTL;
  HAVE_DYNAMIC_BOOST_TEST;
  HAVE_ERT;
  HAVE_SUITESPARSE_UMFPACK_H;
  HAVE_NULLPTR;
  HAVE_STATIC_ASSERT
  ")
include (UseDynamicBoost)
#debug_find_vars ("opm-core")
