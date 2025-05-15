/*
  Copyright 2024 SINTEF AS

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

#ifndef OPM_CUISTL_DEVICE_BLOCK_OPERATIONS_HPP
#define OPM_CUISTL_DEVICE_BLOCK_OPERATIONS_HPP

#include <config.h>
#include <cuda_runtime.h>

#include <cassert>

/*
    This file provides inlineable functions intended for CUDA kernels operating on block matrix elements
    The functions provides various matrix operations that are used by the preconditioners.
*/

namespace
{
// TODO: figure out if this can be generalized effectively, this seems excessively verbose
// explicit formulas based on Dune cpu code
template <class T, int blocksize>
__device__ __forceinline__ void
invBlockOutOfPlace(const T* __restrict__ srcBlock, T* __restrict__ dstBlock)
{
    if (blocksize == 1) {
        dstBlock[0] = T(1.0) / (srcBlock[0]);
    } else if (blocksize == 2) {
        T detInv = T(1.0) / (srcBlock[0] * srcBlock[3] - srcBlock[1] * srcBlock[2]);
        dstBlock[0] = srcBlock[3] * detInv;
        dstBlock[1] = -srcBlock[1] * detInv;
        dstBlock[2] = -srcBlock[2] * detInv;
        dstBlock[3] = srcBlock[0] * detInv;
    } else if (blocksize == 3) {
        // based on Dune implementation
        T t4 = srcBlock[0] * srcBlock[4];
        T t6 = srcBlock[0] * srcBlock[5];
        T t8 = srcBlock[1] * srcBlock[3];
        T t10 = srcBlock[2] * srcBlock[3];
        T t12 = srcBlock[1] * srcBlock[6];
        T t14 = srcBlock[2] * srcBlock[6];

        T t17 = T(1.0)
            / (t4 * srcBlock[8] - t6 * srcBlock[7] - t8 * srcBlock[8] + t10 * srcBlock[7] + t12 * srcBlock[5]
               - t14 * srcBlock[4]); // t17 is 1/determinant

        dstBlock[0] = (srcBlock[4] * srcBlock[8] - srcBlock[5] * srcBlock[7]) * t17;
        dstBlock[1] = -(srcBlock[1] * srcBlock[8] - srcBlock[2] * srcBlock[7]) * t17;
        dstBlock[2] = (srcBlock[1] * srcBlock[5] - srcBlock[2] * srcBlock[4]) * t17;
        dstBlock[3] = -(srcBlock[3] * srcBlock[8] - srcBlock[5] * srcBlock[6]) * t17;
        dstBlock[4] = (srcBlock[0] * srcBlock[8] - t14) * t17;
        dstBlock[5] = -(t6 - t10) * t17;
        dstBlock[6] = (srcBlock[3] * srcBlock[7] - srcBlock[4] * srcBlock[6]) * t17;
        dstBlock[7] = -(srcBlock[0] * srcBlock[7] - t12) * t17;
        dstBlock[8] = (t4 - t8) * t17;
    }
}

// explicit formulas based on Dune cpu code
template <class T, int blocksize>
__device__ __forceinline__ void
invBlockInPlace(T* __restrict__ block)
{
    if (blocksize == 1) {
        block[0] = T(1.0) / (block[0]);
    } else if (blocksize == 2) {
        T detInv = T(1.0) / (block[0] * block[3] - block[1] * block[2]);

        T temp = block[0];
        block[0] = block[3] * detInv;
        block[1] = -block[1] * detInv;
        block[2] = -block[2] * detInv;
        block[3] = temp * detInv;
    } else if (blocksize == 3) {
        T t4 = block[0] * block[4];
        T t6 = block[0] * block[5];
        T t8 = block[1] * block[3];
        T t10 = block[2] * block[3];
        T t12 = block[1] * block[6];
        T t14 = block[2] * block[6];

        T det = (t4 * block[8] - t6 * block[7] - t8 * block[8] + t10 * block[7] + t12 * block[5] - t14 * block[4]);
        T t17 = T(1.0) / det;

        T matrix01 = block[1];
        T matrix00 = block[0];
        T matrix10 = block[3];
        T matrix11 = block[4];

        block[0] = (block[4] * block[8] - block[5] * block[7]) * t17;
        block[1] = -(block[1] * block[8] - block[2] * block[7]) * t17;
        block[2] = (matrix01 * block[5] - block[2] * block[4]) * t17;
        block[3] = -(block[3] * block[8] - block[5] * block[6]) * t17;
        block[4] = (matrix00 * block[8] - t14) * t17;
        block[5] = -(t6 - t10) * t17;
        block[6] = (matrix10 * block[7] - matrix11 * block[6]) * t17;
        block[7] = -(matrix00 * block[7] - t12) * t17;
        block[8] = (t4 - t8) * t17;
    }
}

enum class MVType { SET, PLUS, MINUS };
// SET:   c  = A*b
// PLS:   c += A*b
// MINUS: c -= A*b
template <class T, int blocksize, MVType OP>
__device__ __forceinline__ void
matrixVectorProductWithAction(const T* A, const T* b, T* c)
{
    for (int i = 0; i < blocksize; ++i) {
        if (OP == MVType::SET) {
            c[i] = 0;
        }

        for (int j = 0; j < blocksize; ++j) {
            if (OP == MVType::SET || OP == MVType::PLUS) {
                c[i] += A[i * blocksize + j] * b[j];
            } else if (OP == MVType::MINUS) {
                c[i] -= A[i * blocksize + j] * b[j];
            }
        }
    }
}

template <class T, int blocksize>
__device__ __forceinline__ void
mv(const T* a, const T* b, T* c)
{
    matrixVectorProductWithAction<T, blocksize, MVType::SET>(a, b, c);
}

template <class T, int blocksize>
__device__ __forceinline__ void
umv(const T* a, const T* b, T* c)
{
    matrixVectorProductWithAction<T, blocksize, MVType::PLUS>(a, b, c);
}

template <class T, int blocksize>
__device__ __forceinline__ void
mmv(const T* a, const T* b, T* c)
{
    matrixVectorProductWithAction<T, blocksize, MVType::MINUS>(a, b, c);
}

// dst -= A*B*C
template <class T, int blocksize>
__device__ __forceinline__ void
mmx2Subtraction(const T* A, const T* B, const T* C, T* dst)
{

    T tmp[blocksize * blocksize] = {0};
    // tmp = A*B
    for (int i = 0; i < blocksize; ++i) {
        for (int k = 0; k < blocksize; ++k) {
            for (int j = 0; j < blocksize; ++j) {
                tmp[i * blocksize + j] += A[i * blocksize + k] * B[k * blocksize + j];
            }
        }
    }

    // dst = tmp*C
    for (int i = 0; i < blocksize; ++i) {
        for (int k = 0; k < blocksize; ++k) {
            for (int j = 0; j < blocksize; ++j) {
                dst[i * blocksize + j] -= tmp[i * blocksize + k] * C[k * blocksize + j];
            }
        }
    }
}

// C = A*B, assumes the three buffers do not overlap
template <class T, int blocksize>
__device__ __forceinline__ void
mmNoOverlap(T* A, T* B, T* C)
{
    for (int i = 0; i < blocksize; ++i) {
        for (int k = 0; k < blocksize; ++k) {
            for (int j = 0; j < blocksize; ++j) {
                C[i * blocksize + j] += A[i * blocksize + k] * B[k * blocksize + j];
            }
        }
    }
}

// C = A*B, buffers may overlap
template <class T, int blocksize>
__device__ __forceinline__ void
mmOverlap(T* A, T* B, T* C)
{
    T tmp[blocksize * blocksize] = {0};
    for (int i = 0; i < blocksize; ++i) {
        for (int k = 0; k < blocksize; ++k) {
            for (int j = 0; j < blocksize; ++j) {
                tmp[i * blocksize + j] += A[i * blocksize + k] * B[k * blocksize + j];
            }
        }
    }

    for (int i = 0; i < blocksize; ++i) {
        for (int j = 0; j < blocksize; ++j) {
            C[i * blocksize + j] = tmp[i * blocksize + j];
        }
    }
}

// A -= B
template <class T, int blocksize>
__device__ __forceinline__ void
matrixSubtraction(T* A, T* B)
{

    for (int i = 0; i < blocksize; ++i) {
        for (int j = 0; j < blocksize; ++j) {
            A[i * blocksize + j] -= B[i * blocksize + j];
        }
    }
}

// B = A
template<int blocksize, class ScalarInputType, class ScalarOutputType>
__device__ __forceinline__ void
moveBlock(const ScalarInputType* __restrict__ A, ScalarOutputType* __restrict__ B){
    for (int i = 0; i < blocksize; ++i){
        for (int j = 0; j < blocksize; ++j){
            B[i * blocksize + j] = ScalarOutputType(A[i * blocksize + j]);
        }
    }
}

// TODO: consider merging with existing block operations
// mixed precision general version of c = Ab
template <int blocksize, class MatrixScalar, class VectorScalar, class ResultScalar, class ComputeScalar>
__device__ __forceinline__ void
mvMixedGeneral(const MatrixScalar* A, const VectorScalar* b, ResultScalar* c)
{
    for (int i = 0; i < blocksize; ++i) {
        c[i] = 0;

        for (int j = 0; j < blocksize; ++j) {
            c[i] += ResultScalar(ComputeScalar(A[i * blocksize + j]) * ComputeScalar(b[j]));
        }
    }
}

// TODO: consider merging with existing block operations
// mixed precision general version of c += Ab
template <int blocksize, class MatrixScalar, class VectorScalar, class ResultScalar, class ComputeScalar>
__device__ __forceinline__ void
umvMixedGeneral(const MatrixScalar* A, const VectorScalar* b, ResultScalar* c)
{
    for (int i = 0; i < blocksize; ++i) {
        for (int j = 0; j < blocksize; ++j) {
            c[i] += ResultScalar(ComputeScalar(A[i * blocksize + j]) * ComputeScalar(b[j]));
        }
    }
}

// TODO: consider merging with existing block operations
// Mixed precision general version of c -= Ab
template <int blocksize, class MatrixScalar, class VectorScalar, class ResultScalar, class ComputeScalar>
__device__ __forceinline__ void
mmvMixedGeneral(const MatrixScalar* A, const VectorScalar* b, ResultScalar* c)
{
    for (int i = 0; i < blocksize; ++i) {
        for (int j = 0; j < blocksize; ++j) {
            c[i] -= ResultScalar(ComputeScalar(A[i * blocksize + j]) * ComputeScalar(b[j]));
        }
    }
}

// Checks if a value is close to zero based on a precision limit
// Tolerance is 1e-40, to match matrixblock.hh
template <typename T>
__device__ __forceinline__ bool
isCloseToZero(const T value, const T limit = T(1e-40))
{
    return abs(value) < limit;
}

// Solve a linear system Ax=b for block sizes 1-3
template <typename T, int blockSize>
__device__ __forceinline__ bool
solveBlock(const T* A, const T* b, T* x)
{
    if constexpr (blockSize == 1) {
        if (isCloseToZero(A[0])) {
            return false;
        }
        x[0] = b[0] / A[0];
        return true;
    }
    else if constexpr (blockSize == 2) {
        // Calculate determinant
        T det = A[0] * A[3] - A[1] * A[2];

        if (isCloseToZero(det)) {
            return false;
        }

        T invDet = T(1.0) / det;

        // Compute solution using Cramer's rule
        x[0] = (A[3] * b[0] - A[1] * b[1]) * invDet;
        x[1] = (A[0] * b[1] - A[2] * b[0]) * invDet;

        return true;
    }
    else if constexpr (blockSize == 3) {
        // Calculate determinant
        T det = A[0] * (A[4] * A[8] - A[5] * A[7]) -
                A[1] * (A[3] * A[8] - A[5] * A[6]) +
                A[2] * (A[3] * A[7] - A[4] * A[6]);

        if (isCloseToZero(det)) {
            return false;
        }

        T invDet = T(1.0) / det;

        // Calculate cofactors for each element of b
        x[0] = ((A[4] * A[8] - A[5] * A[7]) * b[0] +
                (A[2] * A[7] - A[1] * A[8]) * b[1] +
                (A[1] * A[5] - A[2] * A[4]) * b[2]) * invDet;

        x[1] = ((A[5] * A[6] - A[3] * A[8]) * b[0] +
                (A[0] * A[8] - A[2] * A[6]) * b[1] +
                (A[2] * A[3] - A[0] * A[5]) * b[2]) * invDet;

        x[2] = ((A[3] * A[7] - A[4] * A[6]) * b[0] +
                (A[1] * A[6] - A[0] * A[7]) * b[1] +
                (A[0] * A[4] - A[1] * A[3]) * b[2]) * invDet;

        return true;
    }
    else {
        // Unsupported block size
        return false;
    }
}

// Transpose a block matrix (row-major)
// Note: This function does NOT support in-place transposition (srcBlock == dstBlock)
// The source and destination blocks must be different memory locations
template <class T, int blockSize>
__device__ __forceinline__ void
transposeBlock(const T* srcBlock, T* dstBlock)
{
    assert(srcBlock != dstBlock && "Source and destination blocks must be different");

    for (int i = 0; i < blockSize; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            dstBlock[j * blockSize + i] = srcBlock[i * blockSize + j];
        }
    }
}
} // namespace

#endif
