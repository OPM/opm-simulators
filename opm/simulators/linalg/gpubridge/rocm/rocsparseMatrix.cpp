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

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/gpubridge/rocm/rocsparseMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/BlockedMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/Matrix.hpp>
#include <opm/simulators/linalg/gpubridge/Misc.hpp>

#include <sstream>
#include <iostream>

namespace Opm::Accelerator {

template<class Scalar>
RocmMatrix<Scalar>::
RocmMatrix(int Nb_, 
           int Mb_,
           int nnzbs_,
           unsigned int block_size_)
    : Nb(Nb_),
      Mb(Mb_),
      nnzbs(nnzbs_),
      block_size(block_size_)
{
    HIP_CHECK(hipMalloc((void**)&nnzValues, sizeof(Scalar) * block_size * block_size * nnzbs));
        
    HIP_CHECK(hipMalloc((void**)&colIndices, sizeof(int) * nnzbs));

    HIP_CHECK(hipMalloc((void**)&rowPointers, sizeof(int) * (Nb + 1)));
}

template <class Scalar>
void RocmMatrix<Scalar>::
upload(Scalar *vals,
       int *cols,
       int *rows,
       hipStream_t stream)
{
    HIP_CHECK(hipMemcpyAsync(nnzValues, vals, sizeof(Scalar) * block_size * block_size * nnzbs, hipMemcpyHostToDevice, stream));
    
    HIP_CHECK(hipMemcpyAsync(colIndices, cols, sizeof(int) * nnzbs, hipMemcpyHostToDevice, stream));
    
    HIP_CHECK(hipMemcpyAsync(rowPointers, rows, sizeof(int) * (Nb + 1), hipMemcpyHostToDevice, stream));
}

template <class Scalar>
void RocmMatrix<Scalar>::
upload(Matrix<Scalar> *matrix,
       hipStream_t stream)
{
    if (block_size != 1) {
        OPM_THROW(std::logic_error, "Error trying to upload a BlockedMatrix to RocmMatrix with different block_size");
    }

    upload(matrix->nnzValues.data(), matrix->colIndices.data(), matrix->rowPointers.data(), stream);
}

template <class Scalar>
void RocmMatrix<Scalar>::
upload(BlockedMatrix<Scalar> *matrix,
       hipStream_t stream)
{
    if (matrix->block_size != block_size) {
        OPM_THROW(std::logic_error, "Error trying to upload a BlockedMatrix to RocmMatrix with different block_size");
    }

    upload(matrix->nnzValues, matrix->colIndices, matrix->rowPointers, stream);
}

template <class Scalar>
RocmVector<Scalar>::RocmVector(int N)
    : size(N)
{
    HIP_CHECK(hipMalloc((void**)&nnzValues, sizeof(Scalar) * N));
}

template <class Scalar>
void RocmVector<Scalar>::
upload(Scalar *vals,
       hipStream_t stream) 
{
    HIP_CHECK(hipMemcpyAsync(nnzValues, vals, sizeof(Scalar) * size, hipMemcpyHostToDevice, stream));    
}

#define INSTANTIATE_TYPE(T)       \
    template class RocmVector<T>; \
    template class RocmMatrix<T>;

INSTANTIATE_TYPE(int)
INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
