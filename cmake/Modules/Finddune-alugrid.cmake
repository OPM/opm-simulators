# - Find DUNE ALUgrid library
#
# Defines the following variables:
#   dune-alugrid_INCLUDE_DIRS    Directory of header files
#   dune-alugrid_LIBRARIES       Directory of shared object files
#   dune-alugrid_DEFINITIONS     Defines that must be set to compile
#   dune-alugrid_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_DUNE_ALUGRID            Binary value to use in config.h

# Copyright (C) 2013 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

include (OpmPackage)
find_opm_package (
  # module name
  "dune-alugrid"

  # dependencies
  # TODO: we should probe for all the HAVE_* values listed below;
  # however, we don't actually use them in our implementation, so
  # we just include them to forward here in case anyone else does
  "CXX11Features REQUIRED;
   dune-grid REQUIRED;
   ZLIB REQUIRED;
   METIS
  "
  # header to search for
  "dune/alugrid/grid.hh"

  # library to search for
  "dunealugrid;alugrid_2d;alugrid_parallel;alugrid_serial"

  # defines to be added to compilations
  ""

  # test program
"#include <dune/alugrid/grid.hh>
int main (void) {
   Dune::ALUGrid</*dim=*/2, /*dimWorld=*/2, Dune::simplex, Dune::nonconforming> grid;
   grid.leafGridView().size(/*codim=*/0);
   return 0;
}
"
  # config variables
  "HAVE_DUNE_ALUGRID
  ")

#debug_find_vars ("dune-grid")

# make version number available in config.h
include (UseDuneVer)
find_dune_version ("dune" "alugrid")
