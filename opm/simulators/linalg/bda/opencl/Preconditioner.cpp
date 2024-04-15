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

namespace Opm
{
namespace Accelerator
{


template <unsigned int block_size>
void Preconditioner<block_size>::setOpencl(std::shared_ptr<cl::Context>& context_, std::shared_ptr<cl::CommandQueue>& queue_) {
    context = context_;
    queue = queue_;
}

template <unsigned int block_size>
std::unique_ptr<Preconditioner<block_size>>
Preconditioner<block_size>::create(Type type, bool opencl_ilu_parallel, int verbosity)
{
    switch (type ) {
    case Type::BILU0:
        return std::make_unique<BILU0<block_size> >(opencl_ilu_parallel, verbosity);
    case Type::CPR:
        return std::make_unique<CPR<block_size> >(opencl_ilu_parallel, verbosity);
    case Type::BISAI:
        return std::make_unique<BISAI<block_size> >(opencl_ilu_parallel, verbosity);
    }

    OPM_THROW(std::logic_error,
              "Invalid preconditioner type " + std::to_string(static_cast<int>(type)));
}

template <unsigned int block_size>
bool Preconditioner<block_size>::
analyze_matrix(BlockedMatrix<double>* mat,
               [[maybe_unused]] BlockedMatrix<double>* jacMat)
{
    return analyze_matrix(mat);
}

template <unsigned int block_size>
bool Preconditioner<block_size>::
create_preconditioner(BlockedMatrix<double>* mat,
                      [[maybe_unused]] BlockedMatrix<double>* jacMat)
{
    return create_preconditioner(mat);
}

#define INSTANTIATE_BDA_FUNCTIONS(n)  \
template std::unique_ptr<Preconditioner<n> > Preconditioner<n>::create(Type, bool, int);         \
template void Preconditioner<n>::setOpencl(std::shared_ptr<cl::Context>&, std::shared_ptr<cl::CommandQueue>&); \
template bool Preconditioner<n>::analyze_matrix(BlockedMatrix<double> *, BlockedMatrix<double> *);                             \
template bool Preconditioner<n>::create_preconditioner(BlockedMatrix<double> *, BlockedMatrix<double> *);

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} //namespace Accelerator
} //namespace Opm

