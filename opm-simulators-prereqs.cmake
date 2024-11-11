# defines that must be present in config.h for our headers
set (opm-simulators_CONFIG_VAR
  HAVE_OPM_GRID
  HAVE_PTHREAD
  HAVE_EWOMS
  HAVE_MPI
  HAVE_PETSC
  COMPILE_GPU_BRIDGE
  HAVE_CUDA
  HAVE_OPENCL
  HAVE_OPENCL_HPP
  HAVE_AMGCL
  HAVE_VEXCL
  HAVE_ROCALUTION
  HAVE_ROCSPARSE
  HAVE_SUITESPARSE_UMFPACK_H
  HAVE_DUNE_COMMON
  HAVE_DUNE_ISTL
  DUNE_ISTL_WITH_CHECKING
  DUNE_ISTL_VERSION_MAJOR
  DUNE_ISTL_VERSION_MINOR
  DUNE_ISTL_VERSION_REVISION
  HAVE_SUITESPARSE_UMFPACK
  HAVE_DAMARIS
  HAVE_HDF5
  USE_HIP
  USE_TRACY
  FLOW_INSTANTIATE_FLOAT
  HAVE_FLOATING_POINT_FROM_CHARS
  OPM_COMPILE_COMPONENTS_TEMPLATE_LIST
  )

# dependencies
set (opm-simulators_DEPS
  # Compile with C99 support if available
  "C99"
  # Various runtime library enhancements
  "Boost 1.44.0
    COMPONENTS date_time system unit_test_framework REQUIRED"
  # DUNE prerequisites
  "dune-common REQUIRED"
  "dune-istl REQUIRED"
  "dune-alugrid"
  "dune-fem"
  # matrix library
  "BLAS REQUIRED"
  "LAPACK REQUIRED"
  # Look for MPI support
  "MPI"
  # Tim Davis' SuiteSparse archive
  "SuiteSparse REQUIRED COMPONENTS umfpack"
  # SuperLU direct solver
  "SuperLU"
  # ROCALUTION from ROCM framework
  "rocalution"
  # packages from ROCm framework
  "rocblas"
  "rocsparse"
  # OPM dependency
  "opm-common REQUIRED"
  "opm-grid REQUIRED"
  "Damaris 1.9"
  "HDF5"
  "Tracy"
  )

find_package_deps(opm-simulators)

if(NOT HAVE_ECL_INPUT OR NOT HAVE_ECL_OUTPUT)
  message(FATAL_ERROR "Eclipse input/output support required in opm-common")
endif()
