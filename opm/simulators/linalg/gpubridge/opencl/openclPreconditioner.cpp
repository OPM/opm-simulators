/*
  Copyright 2021 Equinor ASA

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

#include <opm/common/TimingMacros.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/gpubridge/opencl/openclBILU0.hpp>
#include <opm/simulators/linalg/gpubridge/opencl/openclBISAI.hpp>
#include <opm/simulators/linalg/gpubridge/opencl/openclCPR.hpp>
#include <opm/simulators/linalg/gpubridge/opencl/openclPreconditioner.hpp>

#include <memory>
#include <string>

namespace Opm::Accelerator {

template<class Scalar, unsigned int block_size>
std::unique_ptr<openclPreconditioner<Scalar,block_size>>
openclPreconditioner<Scalar,block_size>::
create(PreconditionerType type,
       int verbosity,
       bool opencl_ilu_parallel)
{
    switch (type ) {
    case PreconditionerType::BILU0:
        return std::make_unique<openclBILU0<Scalar,block_size>>(opencl_ilu_parallel, verbosity);
    case PreconditionerType::CPR:
        return std::make_unique<openclCPR<Scalar,block_size>>(opencl_ilu_parallel, verbosity);
    case PreconditionerType::BISAI:
        return std::make_unique<openclBISAI<Scalar,block_size>>(opencl_ilu_parallel, verbosity);
    }

    OPM_THROW(std::logic_error,
              "Invalid preconditioner type " + std::to_string(static_cast<int>(type)));
}

template<class Scalar, unsigned int block_size>
void openclPreconditioner<Scalar,block_size>::
setOpencl(std::shared_ptr<cl::Context>& context_,
          std::shared_ptr<cl::CommandQueue>& queue_)
{
    context = context_;
    queue = queue_;
}

#define INSTANTIATE_TYPE(T)                   \
    template class openclPreconditioner<T,1>; \
    template class openclPreconditioner<T,2>; \
    template class openclPreconditioner<T,3>; \
    template class openclPreconditioner<T,4>; \
    template class openclPreconditioner<T,5>; \
    template class openclPreconditioner<T,6>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm::Accelerator
