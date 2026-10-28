/*
  Copyright Equinor ASA 2026

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
#ifndef OPM_GPUSYSTEMPRECONDITIONERFACTORY_HEADER_INCLUDED
#define OPM_GPUSYSTEMPRECONDITIONERFACTORY_HEADER_INCLUDED

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSystemTypes.hpp>
#include <opm/simulators/linalg/is_gpu_operator.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

namespace Opm
{

template <class Operator, class Comm, typename>
struct StandardPreconditioners;

namespace gpuistl
{
    template <typename Scalar>
    using GpuSystemSeqOp = Dune::
        MatrixAdapter<GpuSystemMatrix<Scalar>, GpuSystemVector<Scalar>, GpuSystemVector<Scalar>>;
} // namespace gpuistl

// The generic is_gpu_operator test looks for a GpuVector domain, which the
// system operator does not have.
template <typename Scalar>
struct is_gpu_operator<gpuistl::GpuSystemSeqOp<Scalar>> {
    static constexpr bool value = true;
};

template <>
struct StandardPreconditioners<gpuistl::GpuSystemSeqOp<double>,
                               Dune::Amg::SequentialInformation,
                               void> {
    static void add();
};

#if FLOW_INSTANTIATE_FLOAT
template <>
struct StandardPreconditioners<gpuistl::GpuSystemSeqOp<float>,
                               Dune::Amg::SequentialInformation,
                               void> {
    static void add();
};
#endif

} // namespace Opm

#endif // OPM_GPUSYSTEMPRECONDITIONERFACTORY_HEADER_INCLUDED
