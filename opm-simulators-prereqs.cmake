# defines that must be present in config.h for our headers
set (opm-simulators_CONFIG_VAR
  HAVE_MPI
  COMPILE_GPU_BRIDGE
  HAVE_AVX2_EXTENSION
  HAVE_CUDA
  HAVE_OPENCL
  HAVE_OPENCL_HPP
  HAVE_AMGCL
  HAVE_AMGX
  HAVE_VEXCL
  HAVE_ROCALUTION
  HAVE_ROCSPARSE
  HAVE_SUITESPARSE_UMFPACK
  HAVE_DAMARIS
  HAVE_HDF5
  HAVE_HYPRE
  HAVE_DUNE_ISTL
  HAVE_DUNE_COMMON
  HAVE_DUNE_ALUGRID
  HAVE_DUNE_FEM
  USE_HIP
)

include(CheckAVX2)
check_for_avx2()

# dependencies
set (opm-simulators_DEPS
  # Various runtime library enhancements
  "Boost 1.44.0
    COMPONENTS date_time REQUIRED"
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
  "SuiteSparse REQUIRED COMPONENTS UMFPACK"
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
  "fmt"
)

find_package_deps(opm-simulators)
