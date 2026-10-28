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
#include <config.h>
#include <opm/simulators/linalg/gpuistl/GpuSystemBackend.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSystemPreconditionerFactory.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>

#include <functional>
#include <memory>

// The outer FlexibleSolver over GpuSystemSeqOp and its PreconditionerFactory
// are instantiated in the GPU translation units FlexibleSolver_gpu_instantiate.cpp
// and PreconditionerFactory_gpu_instantiate.cpp.  The CPU well sub-solver is
// instantiated with the CPU system preconditioner in system/SystemPreconditioner.cpp.

namespace
{

template <typename Scalar>
void
addGpuSystemCprSeq()
{
    using O = Opm::gpuistl::GpuSystemSeqOp<Scalar>;
    using F = Opm::PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = Opm::gpuistl::GpuSystemVector<Scalar>;
    using P = Opm::PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op,
                     const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      return std::make_shared<Opm::gpuistl::GpuSystemPreconditioner<Scalar>>(
                          op.getmat(), sysWeightCalc, pressureIndex, prm);
                  });
}

} // namespace

template class Opm::SystemPreconditioner<double,
                                         Opm::gpuistl::GpuSystemBackend<double>::ResOperator,
                                         Dune::Amg::SequentialInformation,
                                         Opm::gpuistl::GpuSystemBackend<double>>;

#if FLOW_INSTANTIATE_FLOAT
template class Opm::SystemPreconditioner<float,
                                         Opm::gpuistl::GpuSystemBackend<float>::ResOperator,
                                         Dune::Amg::SequentialInformation,
                                         Opm::gpuistl::GpuSystemBackend<float>>;
#endif

namespace Opm
{

void
StandardPreconditioners<gpuistl::GpuSystemSeqOp<double>, Dune::Amg::SequentialInformation, void>::
    add()
{
    addGpuSystemCprSeq<double>();
}

#if FLOW_INSTANTIATE_FLOAT
void
StandardPreconditioners<gpuistl::GpuSystemSeqOp<float>, Dune::Amg::SequentialInformation, void>::
    add()
{
    addGpuSystemCprSeq<float>();
}
#endif

} // namespace Opm
