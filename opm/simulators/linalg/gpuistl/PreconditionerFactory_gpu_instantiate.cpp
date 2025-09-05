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
#include <dune/istl/operators.hh>
#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#if HAVE_MPI
#include <opm/simulators/linalg/gpuistl/GpuOwnerOverlapCopy.hpp>
#endif

// NOTE: This is very rudimentary, and will be improved once we
// incorporate MPI in the ISTLSolverGPUISTL class.
template class ::Opm::PreconditionerFactory<Dune::MatrixAdapter<::Opm::gpuistl::GpuSparseMatrixWrapper<double>,
                                                                ::Opm::gpuistl::GpuVector<double>,
                                                                ::Opm::gpuistl::GpuVector<double>>,
                                            ::Opm::CommSeq>;

#if FLOW_INSTANTIATE_FLOAT
template class ::Opm::PreconditionerFactory<Dune::MatrixAdapter<::Opm::gpuistl::GpuSparseMatrixWrapper<float>,
                                                                ::Opm::gpuistl::GpuVector<float>,
                                                                ::Opm::gpuistl::GpuVector<float>>,
                                            ::Opm::CommSeq>;
#endif

#if HAVE_MPI
template <class realtype>
using CommGpu = ::Opm::gpuistl::GpuOwnerOverlapCopy<realtype, ::Opm::CommPar>;

template <class Scalar>
using ParOpGpu = Dune::OverlappingSchwarzOperator<::Opm::gpuistl::GpuSparseMatrixWrapper<Scalar>,
                                                  ::Opm::gpuistl::GpuVector<Scalar>,
                                                  ::Opm::gpuistl::GpuVector<Scalar>,
                                                  CommGpu<Scalar>>;

template class ::Opm::PreconditionerFactory<ParOpGpu<double>, CommGpu<double>>;
template class ::Opm::PreconditionerFactory<ParOpGpu<double>, ::Opm::CommSeq>;

#if FLOW_INSTANTIATE_FLOAT
template class ::Opm::PreconditionerFactory<ParOpGpu<float>, CommGpu<float>>;
template class ::Opm::PreconditionerFactory<ParOpGpu<float>, ::Opm::CommSeq>;
#endif
#endif
