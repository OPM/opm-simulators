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
#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fmt/core.h>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cublas_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cublas_wrapper.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_constants.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_pointer_attributes.hpp>
#include <opm/simulators/linalg/gpuistl/detail/vector_operations.hpp>

namespace Opm::gpuistl
{

template <class T>
GpuVector<T>::GpuVector(const std::vector<T>& data)
    : GpuVector(data.data(), data.size())
{
}

template <class T>
GpuVector<T>::GpuVector(const size_t numberOfElements)
    : m_buffer(numberOfElements)
{
}

template <class T>
GpuVector<T>::GpuVector(const T* dataOnHost, const size_t numberOfElements)
    : m_buffer(dataOnHost, numberOfElements)
{
}

template <class T>
GpuVector<T>&
GpuVector<T>::operator=(T scalar)
{
    assertHasElements();
    detail::setVectorValue(data(), dim(), scalar);
    return *this;
}

template <class T>
GpuVector<T>&
GpuVector<T>::operator=(const GpuVector<T>& other)
{
    //TODO-H: Call device-to-device copy


    // Only copy data if both vectors have elements and same size
    if (m_buffer.size() > 0 && other.m_buffer.size() > 0) {
        assertSameSize(other);
        OPM_GPU_SAFE_CALL(cudaMemcpy(data(),
                                      other.data(),
                                      dim() * sizeof(T),
                                      cudaMemcpyDeviceToDevice));
    }
    // If both are zero-sized, assignment is trivial (do nothing)
    return *this;
}

template <class T>
GpuVector<T>::GpuVector(const GpuVector<T>& other)
    : m_buffer(other.m_buffer)
{
}

template <typename T>
const T*
GpuVector<T>::data() const
{
    return m_buffer.data();
}

template <typename T>
T*
GpuVector<T>::data()
{
    return m_buffer.data();
}

template <typename T>
typename GpuVector<T>::size_type
GpuVector<T>::dim() const
{
    return m_buffer.size();
}

template <typename T>
void
GpuVector<T>::resize(size_t new_size)
{
    m_buffer.resize(new_size);
}

template <typename T>
std::vector<T>
GpuVector<T>::asStdVector() const
{
    std::vector<T> temporary(dim());
    copyToHost(temporary);
    return temporary;
}

template <typename T>
void
GpuVector<T>::setZeroAtIndexSet(const GpuVector<int>& indexSet)
{
    detail::setZeroAtIndexSet(data(), indexSet.dim(), indexSet.data());
}

template <typename T>
void
GpuVector<T>::assertSameSize(const GpuVector<T>& x) const
{
    assertSameSize(x.dim());
}

template <typename T>
void
GpuVector<T>::assertSameSize(size_t size) const
{
    if (size != dim()) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Given vector has {}, while we have {}.", size, dim()));
    }
}

template <typename T>
void
GpuVector<T>::assertHasElements() const
{
    if (dim() <= 0) {
        OPM_THROW(std::invalid_argument, "We have 0 elements");
    }
}



template <class T>
GpuVector<T>&
GpuVector<T>::operator*=(const T& scalar)
{
    assertHasElements();
    OPM_CUBLAS_SAFE_CALL(detail::cublasScal(m_cuBlasHandle.get(), detail::to_int(dim()), &scalar, data(), 1));
    return *this;
}

template <class T>
GpuVector<T>&
GpuVector<T>::axpy(T alpha, const GpuVector<T>& y)
{
    assertHasElements();
    assertSameSize(y);
    OPM_CUBLAS_SAFE_CALL(detail::cublasAxpy(m_cuBlasHandle.get(), detail::to_int(dim()), &alpha, y.data(), 1, data(), 1));
    return *this;
}

template <class T>
T
GpuVector<T>::dot(const GpuVector<T>& other) const
{
    assertHasElements();
    assertSameSize(other);
    T result = T(0);
    OPM_CUBLAS_SAFE_CALL(
        detail::cublasDot(m_cuBlasHandle.get(), detail::to_int(dim()), data(), 1, other.data(), 1, &result));
    return result;
}
template <class T>
T
GpuVector<T>::two_norm() const
{
    assertHasElements();
    T result = T(0);
    OPM_CUBLAS_SAFE_CALL(detail::cublasNrm2(m_cuBlasHandle.get(), detail::to_int(dim()), data(), 1, &result));
    return result;
}

template <typename T>
T
GpuVector<T>::dot(const GpuVector<T>& other, const GpuVector<int>& indexSet, GpuVector<T>& buffer) const
{
    return detail::innerProductAtIndices(m_cuBlasHandle.get(), data(), other.data(), buffer.data(), indexSet.dim(), indexSet.data());
}

template <typename T>
T
GpuVector<T>::two_norm(const GpuVector<int>& indexSet, GpuVector<T>& buffer) const
{
    // TODO: [perf] Optimize this to a single call
    return std::sqrt(this->dot(*this, indexSet, buffer));
}

template <typename T>
T
GpuVector<T>::dot(const GpuVector<T>& other, const GpuVector<int>& indexSet) const
{
    GpuVector<T> buffer(indexSet.dim());
    return detail::innerProductAtIndices(m_cuBlasHandle.get(), data(), other.data(), buffer.data(), indexSet.dim(), indexSet.data());
}

template <typename T>
T
GpuVector<T>::two_norm(const GpuVector<int>& indexSet) const
{
    GpuVector<T> buffer(indexSet.dim());
    // TODO: [perf] Optimize this to a single call
    return std::sqrt(this->dot(*this, indexSet, buffer));
}
template <class T>
GpuVector<T>&
GpuVector<T>::operator+=(const GpuVector<T>& other)
{
    assertHasElements();
    assertSameSize(other);
    // TODO: [perf] Make a specialized version of this
    return axpy(1.0, other);
}

template <class T>
GpuVector<T>&
GpuVector<T>::operator-=(const GpuVector<T>& other)
{
    assertHasElements();
    assertSameSize(other);
    // TODO: [perf] Make a specialized version of this
    return axpy(-1.0, other);
}


template <class T>
void
GpuVector<T>::copyFromHost(const T* dataPointer, size_t numberOfElements)
{
    if (numberOfElements > dim()) {
        OPM_THROW(std::runtime_error,
                  fmt::format("Requesting to copy too many elements. Vector has {} elements, while {} was requested.",
                              dim(),
                              numberOfElements));
    }
    OPM_GPU_SAFE_CALL(cudaMemcpy(data(), dataPointer, numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
void
GpuVector<T>::copyFromHostAsync(const T* dataPointer, size_t numberOfElements, cudaStream_t stream)
{
    if (numberOfElements > dim()) {
        OPM_THROW(std::runtime_error,
                  fmt::format("Requesting to copy too many elements. Vector has {} elements, while {} was requested.",
                              dim(),
                              numberOfElements));
    }
    // Asynchronous copy. CUDA runtime will use pinned memory if dataPointer is in a registered region.
    OPM_GPU_SAFE_CALL(cudaMemcpyAsync(data(), dataPointer, numberOfElements * sizeof(T), cudaMemcpyHostToDevice, stream));
}

template <class T>
void
GpuVector<T>::copyToHost(T* dataPointer, size_t numberOfElements) const
{
    // Synchronous version: use default stream and then synchronize.
    copyToHostAsync(dataPointer, numberOfElements, detail::DEFAULT_STREAM);
    OPM_GPU_SAFE_CALL(cudaStreamSynchronize(detail::DEFAULT_STREAM));
}

template <class T>
void
GpuVector<T>::copyToHostAsync(T* dataPointer, size_t numberOfElements, cudaStream_t stream) const
{
    assertSameSize(numberOfElements);
    // Asynchronous copy. CUDA runtime will use pinned memory if dataPointer is in a registered region.
    OPM_GPU_SAFE_CALL(cudaMemcpyAsync(dataPointer, data(), numberOfElements * sizeof(T), cudaMemcpyDeviceToHost, stream));
}

template <class T>
void
GpuVector<T>::copyFromHost(const std::vector<T>& vectorOnHost)
{
    copyFromHost(vectorOnHost.data(), vectorOnHost.size());
}

template <class T>
void
GpuVector<T>::copyFromHostAsync(const std::vector<T>& vectorOnHost, cudaStream_t stream)
{
    copyFromHostAsync(vectorOnHost.data(), vectorOnHost.size(), stream);
}

template <class T>
void
GpuVector<T>::copyToHost(std::vector<T>& vectorOnHost) const
{
    copyToHost(vectorOnHost.data(), vectorOnHost.size());
}

template <class T>
void
GpuVector<T>::copyToHostAsync(std::vector<T>& vectorOnHost, cudaStream_t stream) const
{
    copyToHostAsync(vectorOnHost.data(), vectorOnHost.size(), stream);
}

template <class T>
void
GpuVector<T>::copyFromDeviceToDevice(const GpuVector<T>& other)
{
    assertHasElements();
    assertSameSize(other);

    OPM_GPU_SAFE_CALL(cudaMemcpy(data(),
                                other.data(),
                                dim() * sizeof(T),
                                cudaMemcpyDeviceToDevice));
}

template <typename T>
void
GpuVector<T>::prepareSendBuf(GpuVector<T>& buffer, const GpuVector<int>& indexSet) const
{
    return detail::prepareSendBuf(data(), buffer.data(), indexSet.dim(), indexSet.data());
}
template <typename T>
void
GpuVector<T>::syncFromRecvBuf(GpuVector<T>& buffer, const GpuVector<int>& indexSet)
{
    return detail::syncFromRecvBuf(data(), buffer.data(), indexSet.dim(), indexSet.data());
}

template class GpuVector<double>;
template class GpuVector<float>;
template class GpuVector<int>;

} // namespace Opm::gpuistl
