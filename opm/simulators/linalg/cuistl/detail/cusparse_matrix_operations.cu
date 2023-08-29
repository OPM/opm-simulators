/*
  Copyright 2022-2023 SINTEF AS

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
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <stdexcept>
namespace Opm::cuistl::detail
{

namespace
{
    template <class T> __global__ void cuInvertDiagonalAndFlattenBlocksize1(T* matNonZeroValues, int rowIndices[], int colIndices[], size_t numberOfRows, T* vec)
    {
        const auto thrIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 1;

        if (thrIndex < numberOfRows){
            size_t nnzIdx = rowIndices[thrIndex];
            size_t nnzIdxLim = rowIndices[thrIndex+1]; 

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a zero element
            for (;colIndices[nnzIdx] != thrIndex && nnzIdx <= nnzIdxLim;) { 
                ++nnzIdx;
            }

            // pDiagBlock points to the start of where the diagonal block is stored
            T* pDiagBlock = (matNonZeroValues+(blocksize*blocksize*nnzIdx));
            // pVecBlock points to the start of the block element in the vector where the inverse of the diagonal block element should be stored
            T* pVecBlock = (vec + (blocksize*blocksize*thrIndex));

            pVecBlock[0] =   1.0/(pDiagBlock[0]);
        }
    }

    template <class T> __global__ void cuInvertDiagonalAndFlattenBlocksize2(T* matNonZeroValues, int rowIndices[], int colIndices[], size_t numberOfRows, T* vec)
    {
        const auto thrIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 2;

        if (thrIndex < numberOfRows){
            size_t nnzIdx = rowIndices[thrIndex];
            size_t nnzIdxLim = rowIndices[thrIndex+1];

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a zero element
            for (;colIndices[nnzIdx] != thrIndex && nnzIdx <= nnzIdxLim;) {
                ++nnzIdx;
            }

            // pDiagBlock points to the start of where the diagonal block is stored
            T* pDiagBlock = (matNonZeroValues+(blocksize*blocksize*nnzIdx));
            // pVecBlock points to the start of the block element in the vector where the inverse of the diagonal block element should be stored
            T* pVecBlock = (vec + (blocksize*blocksize*thrIndex));

            // based on Dune implementation
            T det_inv = 1.0/(pDiagBlock[0]*pDiagBlock[3] - pDiagBlock[1]*pDiagBlock[2]);
            pVecBlock[0] =   pDiagBlock[3] * det_inv;
            pVecBlock[1] = - pDiagBlock[1] * det_inv;
            pVecBlock[2] = - pDiagBlock[2] * det_inv;
            pVecBlock[3] =   pDiagBlock[0] * det_inv;
        }
    }
    template <class T> __global__ void cuInvertDiagonalAndFlattenBlocksize3(T* matNonZeroValues, int rowIndices[], int colIndices[], size_t numberOfRows, T* vec)
    {
        const auto thrIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 3;

        if (thrIndex < numberOfRows){
            size_t nnzIdx = rowIndices[thrIndex];
            size_t nnzIdxLim = rowIndices[thrIndex+1];

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a zero element
            for (;colIndices[nnzIdx] != thrIndex && nnzIdx <= nnzIdxLim;) {
                ++nnzIdx;
            }

            // pDiagBlock points to the start of where the diagonal block is stored
            T* pDiagBlock = (matNonZeroValues+(blocksize*blocksize*nnzIdx));
            // pVecBlock points to the start of the block element in the vector where the inverse of the diagonal block element should be stored
            T* pVecBlock = (vec + (blocksize*blocksize*thrIndex));

            // based on Dune implementation
            T t4  = pDiagBlock[0] * pDiagBlock[4];
            T t6  = pDiagBlock[0] * pDiagBlock[5];
            T t8  = pDiagBlock[1] * pDiagBlock[3];
            T t10 = pDiagBlock[2] * pDiagBlock[3];
            T t12 = pDiagBlock[1] * pDiagBlock[6];
            T t14 = pDiagBlock[2] * pDiagBlock[6];

            T t17 = 1.0/(t4*pDiagBlock[8]-t6*pDiagBlock[7]-t8*pDiagBlock[8]+
                    t10*pDiagBlock[7]+t12*pDiagBlock[5]-t14*pDiagBlock[4]); // t17 is 1/determinant

            pVecBlock[0] =  (pDiagBlock[4] * pDiagBlock[8] - pDiagBlock[5] * pDiagBlock[7])*t17;
            pVecBlock[1] = -(pDiagBlock[1] * pDiagBlock[8] - pDiagBlock[2] * pDiagBlock[7])*t17;
            pVecBlock[2] =  (pDiagBlock[1] * pDiagBlock[5] - pDiagBlock[2] * pDiagBlock[4])*t17;
            pVecBlock[3] = -(pDiagBlock[3] * pDiagBlock[8] - pDiagBlock[5] * pDiagBlock[6])*t17;
            pVecBlock[4] =  (pDiagBlock[0] * pDiagBlock[8] - t14) * t17;
            pVecBlock[5] = -(t6-t10) * t17;
            pVecBlock[6] =  (pDiagBlock[3] * pDiagBlock[7] - pDiagBlock[4] * pDiagBlock[6]) * t17;
            pVecBlock[7] = -(pDiagBlock[0] * pDiagBlock[7] - t12) * t17;
            pVecBlock[8] =  (t4-t8) * t17;
        }
    }

    constexpr inline size_t getThreads([[maybe_unused]] size_t numberOfRows)
    {
        return 1024;
    }

    inline size_t getBlocks(size_t numberOfRows)
    {
        const auto threads = getThreads(numberOfRows);
        return (numberOfRows + threads - 1) / threads;
    }
} // namespace

template <class T>
void
invertDiagonalAndFlatten(T* d_mat, int rowIndices[], int colIndices[], size_t numberOfRows, size_t blocksize, T* d_vec)
{
    switch (blocksize) {
    case 1:
        cuInvertDiagonalAndFlattenBlocksize1<<<getBlocks(numberOfRows), getThreads(numberOfRows)>>>(d_mat, rowIndices, colIndices, numberOfRows, d_vec);
        break;
    case 2:
        cuInvertDiagonalAndFlattenBlocksize2<<<getBlocks(numberOfRows), getThreads(numberOfRows)>>>(d_mat, rowIndices, colIndices, numberOfRows, d_vec);
        break;
    case 3:
        cuInvertDiagonalAndFlattenBlocksize3<<<getBlocks(numberOfRows), getThreads(numberOfRows)>>>(d_mat, rowIndices, colIndices, numberOfRows, d_vec);
        break;
    default:
        // TODO: Figure out what is why it did not produce an error or any output in the output stream or the DBG file when I forced this case to execute
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not defined for blocksize>3");
        break;
    }
}

template void invertDiagonalAndFlatten(double*, int*, int*, size_t, size_t, double*);
template void invertDiagonalAndFlatten(float*, int*, int*, size_t, size_t, float*);

} // namespace Opm::cuistl::detail
