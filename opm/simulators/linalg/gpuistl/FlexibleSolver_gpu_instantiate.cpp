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
#include "opm/simulators/linalg/FlexibleSolver.hpp"
#include <dune/istl/operators.hh>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

// NOTE: This is very rudimentary, and will be improved once we
// incorporate MPI in the ISTLSolverGPUISTL class.
template class ::Dune::FlexibleSolver<Dune::MatrixAdapter<::Opm::gpuistl::GpuSparseMatrixWrapper<double>,
                                                          ::Opm::gpuistl::GpuVector<double>,
                                                          ::Opm::gpuistl::GpuVector<double>>>;

#if FLOW_INSTANTIATE_FLOAT
template class ::Dune::FlexibleSolver<Dune::MatrixAdapter<::Opm::gpuistl::GpuSparseMatrixWrapper<float>,
                                                          ::Opm::gpuistl::GpuVector<float>,
                                                          ::Opm::gpuistl::GpuVector<float>>>;
#endif

#if HAVE_MPI
template <class realtype>
using CommGpu = ::Opm::gpuistl::GpuOwnerOverlapCopy<realtype, Comm>;

template <class Scalar>
using ParOpGpu = Dune::OverlappingSchwarzOperator<::Opm::gpuistl::GpuSparseMatrixWrapper<Scalar>,
                                                  ::Opm::gpuistl::GpuVector<Scalar>,
                                                  ::Opm::gpuistl::GpuVector<Scalar>,
                                                  CommGpu<Scalar>>;

#define INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(T)                                                                          \
    template class Dune::FlexibleSolver<ParOpGpu<T>>;                                                                  \
    template Dune::FlexibleSolver<ParOpGpu<T>>::FlexibleSolver(                                                        \
        ParOpGpu<T>& op,                                                                                               \
        const CommGpu<T>& comm,                                                                                        \
        const Opm::PropertyTree& prm,                                                                                  \
        const std::function<typename ParOpGpu<T>::domain_type()>& weightsCalculator,                                   \
        std::size_t pressureIndex);

INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_FLEXIBLESOLVER_GPU_MPI(float)
#endif

#endif
