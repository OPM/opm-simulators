# defines that must be present in config.h for our headers
set (opm-simulators_CONFIG_VAR
  HAVE_OPM_GRID
  HAVE_PTHREAD
  HAVE_EWOMS
  HAVE_ERT
  HAVE_MPI
  HAVE_PETSC
  HAVE_SUITESPARSE_UMFPACK_H
  HAVE_DUNE_ISTL
  DUNE_ISTL_VERSION_MAJOR
  DUNE_ISTL_VERSION_MINOR
  DUNE_ISTL_VERSION_REVISION
  HAVE_SUITESPARSE_UMFPACK
  )

# dependencies
set (opm-simulators_DEPS
  # Compile with C99 support if available
  "C99"
  # Compile with C++0x/11 support if available
  "CXX11Features"
  # Various runtime library enhancements
  "Boost 1.44.0
    COMPONENTS date_time filesystem system unit_test_framework serialization REQUIRED"
  # DUNE prerequisites
  "dune-common REQUIRED"
  "dune-istl REQUIRED"
  # matrix library
  "BLAS REQUIRED"
  "LAPACK REQUIRED"
  # Look for MPI support
  "MPI"
  # Tim Davis' SuiteSparse archive
  "SuiteSparse COMPONENTS umfpack"
  # SuperLU direct solver
  "SuperLU"
  # OPM dependency
  "opm-common REQUIRED"
  "opm-material REQUIRED"
  "opm-grid REQUIRED"
  "ewoms REQUIRED"
  )

find_package_deps(opm-simulators)

if(NOT HAVE_ECL_INPUT OR NOT HAVE_ECL_OUTPUT)
  message(FATAL_ERROR "Eclipse input/output support required in opm-common")
endif()
