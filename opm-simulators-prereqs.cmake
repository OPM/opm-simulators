# defines that must be present in config.h for our headers
set (opm-simulators_CONFIG_VAR
  HAVE_OPM_GRID
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
  HAVE_HYPRE
  USE_HIP
  USE_TRACY
  FLOW_INSTANTIATE_FLOAT
  HAVE_FLOATING_POINT_FROM_CHARS
  OPM_COMPILE_COMPONENTS_TEMPLATE_LIST
)

# CMake 3.30.0 requires to find Boost in CONFIG mode
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30.0)
	set(_Boost_CONFIG_MODE CONFIG)
endif()

include(CheckAVX2)
check_for_avx2()

# dependencies
set (opm-simulators_DEPS
  # Various runtime library enhancements
  "Boost 1.44.0
    COMPONENTS date_time REQUIRED ${_Boost_CONFIG_MODE}"
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
  "rocblas 3.0"
  "rocsparse 4.0"
  # OPM dependency
  "opm-common REQUIRED"
  "opm-grid REQUIRED"
  "Damaris 1.9"
  "HDF5"
  "fmt"
  )

find_package_deps(opm-simulators)

if(NOT HAVE_ECL_INPUT OR NOT HAVE_ECL_OUTPUT)
  message(FATAL_ERROR "Eclipse input/output support required in opm-common")
endif()
