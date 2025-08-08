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
#include <opm/simulators/linalg/gpuistl/detail/vector_operations.hpp>

namespace Opm::gpuistl
{

template <class T>
GpuVector<T>::GpuVector(const std::vector<T>& data)
    : GpuVector(data.data(), detail::to_int(data.size()))
{
}

template <class T>
GpuVector<T>::GpuVector(const size_t numberOfElements)
    : m_numberOfElements(detail::to_int(numberOfElements))
    , m_cuBlasHandle(detail::CuBlasHandle::getInstance())
{
    OPM_GPU_SAFE_CALL(cudaMalloc(&m_dataOnDevice, sizeof(T) * detail::to_size_t(m_numberOfElements)));
}

template <class T>
GpuVector<T>::GpuVector(const T* dataOnHost, const size_t numberOfElements)
    : GpuVector(numberOfElements)
{

    OPM_GPU_SAFE_CALL(cudaMemcpy(
        m_dataOnDevice, dataOnHost, detail::to_size_t(m_numberOfElements) * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
GpuVector<T>&
GpuVector<T>::operator=(T scalar)
{
    assertHasElements();
    detail::setVectorValue(data(), detail::to_size_t(m_numberOfElements), scalar);
    return *this;
}

template <class T>
GpuVector<T>&
GpuVector<T>::operator=(const GpuVector<T>& other)
{
    // Only copy data if both vectors have elements and same size
    if (m_numberOfElements > 0 && other.m_numberOfElements > 0) {
        assertSameSize(other);
        OPM_GPU_SAFE_CALL(cudaMemcpy(m_dataOnDevice,
                                      other.m_dataOnDevice,
                                      detail::to_size_t(m_numberOfElements) * sizeof(T),
                                      cudaMemcpyDeviceToDevice));
    }
    // If both are zero-sized, assignment is trivial (do nothing)
    return *this;
}

template <class T>
GpuVector<T>::GpuVector(const GpuVector<T>& other)
    : GpuVector(other.m_numberOfElements)
{
    // Only copy data if both vectors have elements and same size
    if (m_numberOfElements > 0) {
        assertSameSize(other);
        OPM_GPU_SAFE_CALL(cudaMemcpy(m_dataOnDevice,
                                      other.m_dataOnDevice,
                                      detail::to_size_t(m_numberOfElements) * sizeof(T),
                                      cudaMemcpyDeviceToDevice));
    }
    // If other is zero-sized, assignment is trivial (do nothing)
}

template <class T>
GpuVector<T>::~GpuVector()
{
    OPM_GPU_WARN_IF_ERROR(cudaFree(m_dataOnDevice));
}

template <typename T>
const T*
GpuVector<T>::data() const
{
    return m_dataOnDevice;
}

template <typename T>
typename GpuVector<T>::size_type
GpuVector<T>::dim() const
{
    // Note that there is no way for m_numberOfElements to be non-positive,
    // but for sanity we still use the safe conversion function here.
    //
    // We also doubt that this will lead to any performance penality, but should this prove
    // to be false, this can be replaced by a simple cast to size_t
    return detail::to_size_t(m_numberOfElements);
}

template <typename T>
void
GpuVector<T>::resize(size_t new_size)
{
    const int new_elements = detail::to_int(new_size);

    if (new_elements == m_numberOfElements) {
        return;
    }

    if (new_elements == 0) {
        // Free existing memory and set to empty state
        if (m_dataOnDevice != nullptr) {
            OPM_GPU_WARN_IF_ERROR(cudaFree(m_dataOnDevice));
            m_dataOnDevice = nullptr;
        }
        m_numberOfElements = 0;
        return;
    }

    // Allocate new memory
    T* new_data = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&new_data, sizeof(T) * new_size));

    if (m_dataOnDevice != nullptr && m_numberOfElements > 0) {
        // Copy existing data (up to the minimum of old and new size)
        const size_t copy_elements = std::min(detail::to_size_t(m_numberOfElements), new_size);
        if (copy_elements > 0) {
            OPM_GPU_SAFE_CALL(cudaMemcpy(new_data, m_dataOnDevice,
                                          sizeof(T) * copy_elements,
                                          cudaMemcpyDeviceToDevice));
        }

        // Free old memory
        OPM_GPU_WARN_IF_ERROR(cudaFree(m_dataOnDevice));
    }

    m_dataOnDevice = new_data;
    m_numberOfElements = new_elements;
}

template <typename T>
std::vector<T>
GpuVector<T>::asStdVector() const
{
    std::vector<T> temporary(detail::to_size_t(m_numberOfElements));
    copyToHost(temporary);
    return temporary;
}

template <typename T>
void
GpuVector<T>::setZeroAtIndexSet(const GpuVector<int>& indexSet)
{
    detail::setZeroAtIndexSet(m_dataOnDevice, indexSet.dim(), indexSet.data());
}

template <typename T>
void
GpuVector<T>::assertSameSize(const GpuVector<T>& x) const
{
    assertSameSize(x.m_numberOfElements);
}

template <typename T>
void
GpuVector<T>::assertSameSize(int size) const
{
    if (size != m_numberOfElements) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Given vector has {}, while we have {}.", size, m_numberOfElements));
    }
}

template <typename T>
void
GpuVector<T>::assertHasElements() const
{
    if (m_numberOfElements <= 0) {
        OPM_THROW(std::invalid_argument, "We have 0 elements");
    }
}

template <typename T>
T*
GpuVector<T>::data()
{
    return m_dataOnDevice;
}

template <class T>
GpuVector<T>&
GpuVector<T>::operator*=(const T& scalar)
{
    assertHasElements();
    OPM_CUBLAS_SAFE_CALL(detail::cublasScal(m_cuBlasHandle.get(), m_numberOfElements, &scalar, data(), 1));
    return *this;
}

template <class T>
GpuVector<T>&
GpuVector<T>::axpy(T alpha, const GpuVector<T>& y)
{
    assertHasElements();
    assertSameSize(y);
    OPM_CUBLAS_SAFE_CALL(detail::cublasAxpy(m_cuBlasHandle.get(), m_numberOfElements, &alpha, y.data(), 1, data(), 1));
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
        detail::cublasDot(m_cuBlasHandle.get(), m_numberOfElements, data(), 1, other.data(), 1, &result));
    return result;
}
template <class T>
T
GpuVector<T>::two_norm() const
{
    assertHasElements();
    T result = T(0);
    OPM_CUBLAS_SAFE_CALL(detail::cublasNrm2(m_cuBlasHandle.get(), m_numberOfElements, data(), 1, &result));
    return result;
}

template <typename T>
T
GpuVector<T>::dot(const GpuVector<T>& other, const GpuVector<int>& indexSet, GpuVector<T>& buffer) const
{
    return detail::innerProductAtIndices(m_cuBlasHandle.get(), m_dataOnDevice, other.data(), buffer.data(), indexSet.dim(), indexSet.data());
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
    return detail::innerProductAtIndices(m_cuBlasHandle.get(), m_dataOnDevice, other.data(), buffer.data(), indexSet.dim(), indexSet.data());
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
    assertSameSize(detail::to_int(numberOfElements));
    // Asynchronous copy. CUDA runtime will use pinned memory if dataPointer is in a registered region.
    OPM_GPU_SAFE_CALL(cudaMemcpyAsync(dataPointer, data(), numberOfElements * sizeof(T), cudaMemcpyDeviceToHost, stream));
}

template <class T>
void
GpuVector<T>::copyFromHost(const std::vector<T>& data)
{
    copyFromHost(data.data(), data.size());
}

template <class T>
void
GpuVector<T>::copyFromHostAsync(const std::vector<T>& data, cudaStream_t stream)
{
    copyFromHostAsync(data.data(), data.size(), stream);
}

template <class T>
void
GpuVector<T>::copyToHost(std::vector<T>& data) const
{
    copyToHost(data.data(), data.size());
}

template <class T>
void
GpuVector<T>::copyToHostAsync(std::vector<T>& data, cudaStream_t stream) const
{
    copyToHostAsync(data.data(), data.size(), stream);
}

template <class T>
void
GpuVector<T>::copyFromDeviceToDevice(const GpuVector<T>& data) const
{
    assertHasElements();
    assertSameSize(data);

    OPM_GPU_SAFE_CALL(cudaMemcpy(m_dataOnDevice,
                                data.m_dataOnDevice,
                                detail::to_size_t(m_numberOfElements) * sizeof(T),
                                cudaMemcpyDeviceToDevice));
}

template <typename T>
void
GpuVector<T>::prepareSendBuf(GpuVector<T>& buffer, const GpuVector<int>& indexSet) const
{
    return detail::prepareSendBuf(m_dataOnDevice, buffer.data(), indexSet.dim(), indexSet.data());
}
template <typename T>
void
GpuVector<T>::syncFromRecvBuf(GpuVector<T>& buffer, const GpuVector<int>& indexSet) const
{
    return detail::syncFromRecvBuf(m_dataOnDevice, buffer.data(), indexSet.dim(), indexSet.data());
}

template class GpuVector<double>;
template class GpuVector<float>;
template class GpuVector<int>;

} // namespace Opm::gpuistl
