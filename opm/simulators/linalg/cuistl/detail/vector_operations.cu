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
#include <opm/simulators/linalg/cuistl/detail/vector_operations.hpp>
// TODO: [perf] Get rid of thrust.
#include <stdexcept>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
namespace Opm::cuistl::detail
{

namespace
{
    template <class T>
    __global__ void setZeroAtIndexSetKernel(T* devicePointer, size_t numberOfElements, const int* indices)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;

        if (globalIndex < numberOfElements) {
            devicePointer[indices[globalIndex]] = T(0);
        }
    }

    template <class T>
    __global__ void setVectorValueKernel(T* devicePointer, size_t numberOfElements, T value)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;

        if (globalIndex < numberOfElements) {
            devicePointer[globalIndex] = value;
        }
    }

    template <class T>
    __global__ void elementWiseMultiplyMVKernelBlocksize1(T* squareBlockVector, const size_t numberOfElements, T* vec)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 1;

        if (globalIndex < numberOfElements) {
            T* pMat = (squareBlockVector + (blocksize * blocksize * globalIndex));
            T* pVec = (vec + (blocksize * globalIndex));

            T v0 = pVec[0];
            pVec[0] = pMat[0] * v0;
        }
    }

    template <class T>
    __global__ void elementWiseMultiplyMVKernelBlocksize2(T* squareBlockVector, const size_t numberOfElements, T* vec)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 2;

        if (globalIndex < numberOfElements) {
            T* pMat = (squareBlockVector + (blocksize * blocksize * globalIndex));
            T* pVec = (vec + (blocksize * globalIndex));

            T v0 = pVec[0], v1 = pVec[1];
            pVec[0] = pMat[0] * v0 + pMat[1] * v1;
            pVec[1] = pMat[2] * v0 + pMat[3] * v1;
        }
    }

    template <class T>
    __global__ void elementWiseMultiplyMVKernelBlocksize3(T* squareBlockVector, const size_t numberOfElements, T* vec)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 3;

        if (globalIndex < numberOfElements) {
            T* pMat = (squareBlockVector + (blocksize * blocksize * globalIndex));
            T* pVec = (vec + (blocksize * globalIndex));

            T v0 = pVec[0], v1 = pVec[1], v2 = pVec[2];
            pVec[0] = pMat[0] * v0 + pMat[1] * v1 + pMat[2] * v2;
            pVec[1] = pMat[3] * v0 + pMat[4] * v1 + pMat[5] * v2;
            pVec[2] = pMat[6] * v0 + pMat[7] * v1 + pMat[8] * v2;
        }
    }

    template <class T>
    __global__ void
    elementWiseMultiplyKernel(const T* a, const T* b, T* buffer, size_t numberOfElements, const int* indices)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;

        // TODO: [perf] Is it faster to just use a mask? Probably does not matter either way
        //              This is hopefully not where we will spend most of our time.
        if (globalIndex < numberOfElements) {
            buffer[globalIndex] = a[indices[globalIndex]] * b[indices[globalIndex]];
        }
    }

    constexpr inline size_t getThreads([[maybe_unused]] size_t numberOfElements)
    {
        return 1024;
    }

    inline size_t getBlocks(size_t numberOfElements)
    {
        const auto threads = getThreads(numberOfElements);
        return (numberOfElements + threads - 1) / threads;
    }
} // namespace

template <class T>
void
setVectorValue(T* deviceData, size_t numberOfElements, const T& value)
{
    setVectorValueKernel<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
        deviceData, numberOfElements, value);
}

template void setVectorValue(double*, size_t, const double&);
template void setVectorValue(float*, size_t, const float&);
template void setVectorValue(int*, size_t, const int&);

template <class T>
void
setZeroAtIndexSet(T* deviceData, size_t numberOfElements, const int* indices)
{
    setZeroAtIndexSetKernel<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
        deviceData, numberOfElements, indices);
}
template void setZeroAtIndexSet(double*, size_t, const int*);
template void setZeroAtIndexSet(float*, size_t, const int*);
template void setZeroAtIndexSet(int*, size_t, const int*);

template <class T>
T
innerProductAtIndices(const T* deviceA, const T* deviceB, T* buffer, size_t numberOfElements, const int* indices)
{
    elementWiseMultiplyKernel<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
        deviceA, deviceB, buffer, numberOfElements, indices);

    // TODO: [perf] Get rid of thrust and use a more direct reduction here.
    auto bufferAsDevicePointer = thrust::device_pointer_cast(buffer);
    return thrust::reduce(bufferAsDevicePointer, bufferAsDevicePointer + numberOfElements);
}

template double innerProductAtIndices(const double*, const double*, double* buffer, size_t, const int*);
template float innerProductAtIndices(const float*, const float*, float* buffer, size_t, const int*);
template int innerProductAtIndices(const int*, const int*, int* buffer, size_t, const int*);


template <class T>
void
blockVectorMultiplicationAtAllIndices(T* squareBlockVector,
                                      const size_t numberOfElements,
                                      const size_t blocksize,
                                      T* vec)
{
    switch (blocksize) {
    case 1:
        elementWiseMultiplyMVKernelBlocksize1<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
            squareBlockVector, numberOfElements, vec);
        break;
    case 2:
        elementWiseMultiplyMVKernelBlocksize2<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
            squareBlockVector, numberOfElements, vec);
        break;
    case 3:
        elementWiseMultiplyMVKernelBlocksize3<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
            squareBlockVector, numberOfElements, vec);
        break;
    default:
        OPM_THROW(std::invalid_argument, "blockvector hadamard product not defined for blocksize>3");
        break;
    }
}

template void blockVectorMultiplicationAtAllIndices(double*, const size_t, const size_t, double*);
template void blockVectorMultiplicationAtAllIndices(float*, const size_t, const size_t, float*);

} // namespace Opm::cuistl::detail
