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

#include <opm/simulators/linalg/gpubridge/rocm/rocsparseBILU0.hpp>
#include <opm/simulators/linalg/gpubridge/rocm/rocsparseCPR.hpp>
#include <opm/simulators/linalg/gpubridge/rocm/rocsparsePreconditioner.hpp>

namespace Opm::Accelerator {

template <class Scalar, unsigned int block_size>
std::unique_ptr<rocsparsePreconditioner<Scalar,block_size> > rocsparsePreconditioner<Scalar,block_size>::
create(PreconditionerType type,
       int verbosity)
{
    switch (type ) {
    case PreconditionerType::BILU0:
        return std::make_unique<Opm::Accelerator::rocsparseBILU0<Scalar, block_size> >(verbosity);
    case PreconditionerType::CPR:
        return std::make_unique<Opm::Accelerator::rocsparseCPR<Scalar, block_size> >(verbosity);
    default:
        OPM_THROW(std::logic_error,
                  "Invalid preconditioner type " + std::to_string(static_cast<int>(type)));
    }
}

template <class Scalar, unsigned int block_size>
void rocsparsePreconditioner<Scalar, block_size>::
set_matrix_analysis(rocsparse_mat_descr desc_L,
                    rocsparse_mat_descr desc_U)
{
    descr_L = desc_L;
    descr_U = desc_U;
}

template <class Scalar, unsigned int block_size>
void rocsparsePreconditioner<Scalar, block_size>::
set_context(rocsparse_handle handle_,
            rocblas_handle blas_handle_,
            rocsparse_direction dir_,
            rocsparse_operation operation_,
            hipStream_t stream_)
{
    this->handle = handle_;
    this->blas_handle = blas_handle_;
    this->dir = dir_;
    this->operation = operation_;
    this->stream = stream_;
}

template <class Scalar, unsigned int block_size>
void rocsparsePreconditioner<Scalar, block_size>::
setJacMat(const BlockedMatrix<Scalar>& jMat)
{
    this->jacMat = std::make_shared<BlockedMatrix<Scalar>>(jMat);
}



#define INSTANTIATE_TYPE(T)                      \
    template class rocsparsePreconditioner<T,1>; \
    template class rocsparsePreconditioner<T,2>; \
    template class rocsparsePreconditioner<T,3>; \
    template class rocsparsePreconditioner<T,4>; \
    template class rocsparsePreconditioner<T,5>; \
    template class rocsparsePreconditioner<T,6>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} //namespace Opm
