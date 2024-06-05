/*
  Copyright 2024 Equinor ASA

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
#include <memory>
#include <opm/common/TimingMacros.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/rocm/rocsparseBILU0.hpp>
#include <opm/simulators/linalg/bda/rocm/rocsparseCPR.hpp>
#include <opm/simulators/linalg/bda/rocm/rocsparsePreconditioner.hpp>

namespace Opm::Accelerator {

template <class Scalar, unsigned int block_size>
std::unique_ptr<rocsparsePreconditioner<Scalar,block_size> > rocsparsePreconditioner<Scalar,block_size>::
create(PreconditionerType type,
       int verbosity)
{
    if (type == PreconditionerType::BILU0) {
        return std::make_unique<Opm::Accelerator::rocsparseBILU0<Scalar, block_size> >(verbosity);
    } else if (type == PreconditionerType::CPR) {
        return std::make_unique<Opm::Accelerator::rocsparseCPR<Scalar, block_size> >(verbosity);
    } else {
        OPM_THROW(std::logic_error, "Invalid PreconditionerType");
    }
}

template <class Scalar, unsigned int block_size>
void rocsparsePreconditioner<Scalar, block_size>::
set_matrix_analysis(rocsparse_mat_descr descr_L,
                    rocsparse_mat_descr descr_U)
{
    descr_L = descr_L;
    descr_U = descr_U;
}

template <class Scalar, unsigned int block_size>
void rocsparsePreconditioner<Scalar, block_size>::
set_context(rocsparse_handle handle,
            rocsparse_direction dir,
            rocsparse_operation operation,
            hipStream_t stream)
{
    this->handle = handle;
    this->dir = dir;
    this->operation = operation;
    this->stream = stream;
}

template <class Scalar, unsigned int block_size>
void rocsparsePreconditioner<Scalar, block_size>::
setJacMat(BlockedMatrix<Scalar> jacMat) {
    this->jacMat = std::make_shared<BlockedMatrix<Scalar>>(jacMat);
}

#define INSTANTIATE_TYPE(T)                  \
    template class rocsparsePreconditioner<T,1>; \
    template class rocsparsePreconditioner<T,2>; \
    template class rocsparsePreconditioner<T,3>; \
    template class rocsparsePreconditioner<T,4>; \
    template class rocsparsePreconditioner<T,5>; \
    template class rocsparsePreconditioner<T,6>;

INSTANTIATE_TYPE(double)

} //namespace Opm

