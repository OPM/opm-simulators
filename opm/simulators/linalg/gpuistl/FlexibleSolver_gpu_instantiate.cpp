/*
  Copyright 2025 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"
#include "dune/istl/schwarz.hh"
#include <dune/istl/operators.hh>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

// NOTE: This is very rudimentary, and will be improved once we
// incorporate MPI in the ISTLSolverGPUISTL class.
template class ::Dune::FlexibleSolver<Dune::MatrixAdapter<::Opm::gpuistl::GpuSparseMatrix<double>,
                                                          ::Opm::gpuistl::GpuVector<double>,
                                                          ::Opm::gpuistl::GpuVector<double>>>;

#if FLOW_INSTANTIATE_FLOAT
template class ::Dune::FlexibleSolver<Dune::MatrixAdapter<::Opm::gpuistl::GpuSparseMatrix<float>,
                                                          ::Opm::gpuistl::GpuVector<float>,
                                                          ::Opm::gpuistl::GpuVector<float>>>;
#endif

#if HAVE_MPI

#define INSTANTIATE_FLEXIBLESOLVER_OP_GPU(COMMTYPE, ...)                                            \
    template class Dune::FlexibleSolver<__VA_ARGS__>;                                               \
    template Dune::FlexibleSolver<__VA_ARGS__>::                                                    \
        FlexibleSolver(__VA_ARGS__& op,                                                             \
                       const COMMTYPE& comm,                                                            \
                       const Opm::PropertyTree& prm,                                                \
                       const std::function<typename __VA_ARGS__::domain_type()>& weightsCalculator, \
                       std::size_t pressureIndex);

#define INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(realtype, blocksize) \
    INSTANTIATE_FLEXIBLESOLVER_OP_GPU(::Opm::gpuistl::GpuOwnerOverlapCopy<realtype, blocksize, Comm>, \
        Dune::OverlappingSchwarzOperator<::Opm::gpuistl::GpuSparseMatrix<realtype>, \
          ::Opm::gpuistl::GpuVector<realtype>, \
          ::Opm::gpuistl::GpuVector<realtype>, \
          ::Opm::gpuistl::GpuOwnerOverlapCopy<realtype, blocksize, Comm>>>);


INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(double, 1)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(double, 2)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(double, 3)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(double, 4)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(double, 5)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(double, 6)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(float, 1)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(float, 2)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(float, 3)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(float, 4)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(float, 5)
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(float, 6)
#endif
#endif
