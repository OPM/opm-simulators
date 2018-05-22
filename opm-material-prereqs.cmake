# defines that must be present in config.h for our headers
set (opm-material_CONFIG_VAR
  HAVE_MPI
  HAVE_TYPE_TRAITS
  HAVE_VALGRIND
  HAVE_FINAL
  HAVE_ECL_INPUT
  )

# dependencies
set (opm-material_DEPS
  # compile with C99 support if available
  "C99"
  # compile with C++0x/11 support if available
  "CXX11Features REQUIRED"
  # prerequisite OPM modules
  "ecl"
  "opm-common REQUIRED"
  # DUNE dependency
  "dune-common REQUIRED"
  # valgrind client requests
  "Valgrind"
  "Boost COMPONENTS filesystem"
  )

find_package_deps(opm-material)
