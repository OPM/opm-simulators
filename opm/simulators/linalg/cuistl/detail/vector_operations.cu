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
#include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cublas_wrapper.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <stdexcept>
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

    template <class T, int blocksize>
    __global__ void weightedDiagMV(
        const T* squareBlockVector, const size_t numberOfElements, const T w, const T* globalSrcVec, T* globalDstVec)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;

        if (globalIndex < numberOfElements) {
            const T* localBlock = (squareBlockVector + (blocksize * blocksize * globalIndex));
            const T* localSrcVec = (globalSrcVec + (blocksize * globalIndex));
            T* localDstVec = (globalDstVec + (blocksize * globalIndex));

            for (int i = 0; i < blocksize; ++i) {
                T rowResult = 0.0;
                for (int j = 0; j < blocksize; ++j) {
                    rowResult += localBlock[i * blocksize + j] * localSrcVec[j];
                }
                localDstVec[i] = rowResult * w;
            }
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

    template <class T>
    __global__ void
    prepareSendBufKernel(const T* a, T* buffer, size_t numberOfElements, const int* indices)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;

        if (globalIndex < numberOfElements) {
            buffer[globalIndex] = a[indices[globalIndex]];
        }
    }
    template <class T>
    __global__ void
    syncFromRecvBufKernel(T* a, T* buffer, size_t numberOfElements, const int* indices)
    {
        const auto globalIndex = blockDim.x * blockIdx.x + threadIdx.x;

        if (globalIndex < numberOfElements) {
            a[indices[globalIndex]] = buffer[globalIndex];
        }
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
innerProductAtIndices(cublasHandle_t cublasHandle, const T* deviceA, const T* deviceB, T* buffer, size_t numberOfElements, const int* indices)
{
    elementWiseMultiplyKernel<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
        deviceA, deviceB, buffer, numberOfElements, indices);

    // TODO: [perf] Get rid of the allocation here.
    CuVector<T> oneVector(numberOfElements);
    oneVector = 1.0;
    T result = 0.0;
    OPM_CUBLAS_SAFE_CALL(cublasDot(cublasHandle, numberOfElements, oneVector.data(), 1, buffer, 1, &result));
    return result;
}

template double innerProductAtIndices(cublasHandle_t, const double*, const double*, double* buffer, size_t, const int*);
template float innerProductAtIndices(cublasHandle_t, const float*, const float*, float* buffer, size_t, const int*);
template int innerProductAtIndices(cublasHandle_t, const int*, const int*, int* buffer, size_t, const int*);

template <class T>
void prepareSendBuf(const T* deviceA, T* buffer, size_t numberOfElements, const int* indices)
{
    prepareSendBufKernel<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(deviceA, buffer, numberOfElements, indices);
    cudaDeviceSynchronize(); // The buffers are prepared for MPI. Wait for them to finish.
}
template void prepareSendBuf(const double* deviceA, double* buffer, size_t numberOfElements, const int* indices);
template void prepareSendBuf(const float* deviceA, float* buffer, size_t numberOfElements, const int* indices);
template void prepareSendBuf(const int* deviceA, int* buffer, size_t numberOfElements, const int* indices);

template <class T>
void syncFromRecvBuf(T* deviceA, T* buffer, size_t numberOfElements, const int* indices)
{
    syncFromRecvBufKernel<<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(deviceA, buffer, numberOfElements, indices);
    //cudaDeviceSynchronize(); // Not needed, I guess...
}
template void syncFromRecvBuf(double* deviceA, double* buffer, size_t numberOfElements, const int* indices);
template void syncFromRecvBuf(float* deviceA, float* buffer, size_t numberOfElements, const int* indices);
template void syncFromRecvBuf(int* deviceA, int* buffer, size_t numberOfElements, const int* indices);

template <class T>
void
weightedDiagMV(const T* squareBlockVector,
               const size_t numberOfElements,
               const size_t blocksize,
               T relaxationFactor,
               const T* srcVec,
               T* dstVec)
{
    switch (blocksize) {
    case 1:
        weightedDiagMV<T, 1><<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
            squareBlockVector, numberOfElements, relaxationFactor, srcVec, dstVec);
        break;
    case 2:
        weightedDiagMV<T, 2><<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
            squareBlockVector, numberOfElements, relaxationFactor, srcVec, dstVec);
        break;
    case 3:
        weightedDiagMV<T, 3><<<getBlocks(numberOfElements), getThreads(numberOfElements)>>>(
            squareBlockVector, numberOfElements, relaxationFactor, srcVec, dstVec);
        break;
    default:
        OPM_THROW(std::invalid_argument, "blockvector Hadamard product not implemented for blocksize>3");
        break;
    }
}

template void weightedDiagMV(const double*, const size_t, const size_t, double, const double*, double*);
template void weightedDiagMV(const float*, const size_t, const size_t, float, const float*, float*);

} // namespace Opm::cuistl::detail
