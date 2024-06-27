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
#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>

#include <opm/common/TimingMacros.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/opencl/BILU0.hpp>
#include <opm/simulators/linalg/bda/opencl/BISAI.hpp>
#include <opm/simulators/linalg/bda/opencl/CPR.hpp>

#include <memory>
#include <string>

namespace Opm::Accelerator {

template<class Scalar, unsigned int block_size>
void Preconditioner<Scalar,block_size>::
 setOpencl(std::shared_ptr<cl::Context>& context_,
           std::shared_ptr<cl::CommandQueue>& queue_)
{
    context = context_;
    queue = queue_;
}

template<class Scalar, unsigned int block_size>
std::unique_ptr<Preconditioner<Scalar,block_size>>
Preconditioner<Scalar,block_size>::create(Type type, bool opencl_ilu_parallel, int verbosity)
{
    switch (type ) {
    case Type::BILU0:
        return std::make_unique<BILU0<Scalar,block_size>>(opencl_ilu_parallel, verbosity);
    case Type::CPR:
        return std::make_unique<CPR<Scalar,block_size>>(opencl_ilu_parallel, verbosity);
    case Type::BISAI:
        return std::make_unique<BISAI<Scalar,block_size>>(opencl_ilu_parallel, verbosity);
    }

    OPM_THROW(std::logic_error,
              "Invalid preconditioner type " + std::to_string(static_cast<int>(type)));
}

template<class Scalar, unsigned int block_size>
bool Preconditioner<Scalar,block_size>::
analyze_matrix(BlockedMatrix<Scalar>* mat,
               [[maybe_unused]] BlockedMatrix<Scalar>* jacMat)
{
    return analyze_matrix(mat);
}

template<class Scalar, unsigned int block_size>
bool Preconditioner<Scalar,block_size>::
create_preconditioner(BlockedMatrix<Scalar>* mat,
                      [[maybe_unused]] BlockedMatrix<Scalar>* jacMat)
{
    return create_preconditioner(mat);
}

#define INSTANTIATE_TYPE(T)             \
    template class Preconditioner<T,1>; \
    template class Preconditioner<T,2>; \
    template class Preconditioner<T,3>; \
    template class Preconditioner<T,4>; \
    template class Preconditioner<T,5>; \
    template class Preconditioner<T,6>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm::Accelerator
