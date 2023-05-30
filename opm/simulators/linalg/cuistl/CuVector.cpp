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
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cublas_wrapper.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/vector_operations.hpp>

namespace Opm::cuistl
{

template <class T>
CuVector<T>::CuVector(const std::vector<T>& data)
    : CuVector(data.data(), detail::to_int(data.size()))
{
}

template <class T>
CuVector<T>::CuVector(const size_t numberOfElements)
    : m_numberOfElements(detail::to_int(numberOfElements))
    , m_cuBlasHandle(detail::CuBlasHandle::getInstance())
{
    OPM_CUDA_SAFE_CALL(cudaMalloc(&m_dataOnDevice, sizeof(T) * detail::to_size_t(m_numberOfElements)));
}

template <class T>
CuVector<T>::CuVector(const T* dataOnHost, const size_t numberOfElements)
    : CuVector(numberOfElements)
{

    OPM_CUDA_SAFE_CALL(cudaMemcpy(
        m_dataOnDevice, dataOnHost, detail::to_size_t(m_numberOfElements) * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
CuVector<T>&
CuVector<T>::operator=(T scalar)
{
    assertHasElements();
    detail::setVectorValue(data(), detail::to_size_t(m_numberOfElements), scalar);
    return *this;
}

template <class T>
CuVector<T>&
CuVector<T>::operator=(const CuVector<T>& other)
{
    assertHasElements();
    assertSameSize(other);

    OPM_CUDA_SAFE_CALL(cudaMemcpy(m_dataOnDevice,
                                  other.m_dataOnDevice,
                                  detail::to_size_t(m_numberOfElements) * sizeof(T),
                                  cudaMemcpyDeviceToDevice));
    return *this;
}

template <class T>
CuVector<T>::CuVector(const CuVector<T>& other)
    : CuVector(other.m_numberOfElements)
{
    assertHasElements();
    assertSameSize(other);
    OPM_CUDA_SAFE_CALL(cudaMemcpy(m_dataOnDevice,
                                  other.m_dataOnDevice,
                                  detail::to_size_t(m_numberOfElements) * sizeof(T),
                                  cudaMemcpyDeviceToDevice));
}

template <class T>
CuVector<T>::~CuVector()
{
    OPM_CUDA_WARN_IF_ERROR(cudaFree(m_dataOnDevice));
}

template <typename T>
const T*
CuVector<T>::data() const
{
    return m_dataOnDevice;
}

template <typename T>
typename CuVector<T>::size_type
CuVector<T>::dim() const
{
    // Note that there is no way for m_numberOfElements to be non-positive,
    // but for sanity we still use the safe conversion function here.
    //
    // We also doubt that this will lead to any performance penality, but should this prove
    // to be false, this can be replaced by a simple cast to size_t
    return detail::to_size_t(m_numberOfElements);
}

template <typename T>
std::vector<T>
CuVector<T>::asStdVector() const
{
    std::vector<T> temporary(detail::to_size_t(m_numberOfElements));
    copyToHost(temporary);
    return temporary;
}

template <typename T>
void
CuVector<T>::setZeroAtIndexSet(const CuVector<int>& indexSet)
{
    detail::setZeroAtIndexSet(m_dataOnDevice, indexSet.dim(), indexSet.data());
}

template <typename T>
void
CuVector<T>::assertSameSize(const CuVector<T>& x) const
{
    assertSameSize(x.m_numberOfElements);
}

template <typename T>
void
CuVector<T>::assertSameSize(int size) const
{
    if (size != m_numberOfElements) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Given vector has {}, while we have {}.", size, m_numberOfElements));
    }
}

template <typename T>
void
CuVector<T>::assertHasElements() const
{
    if (m_numberOfElements <= 0) {
        OPM_THROW(std::invalid_argument, "We have 0 elements");
    }
}

template <typename T>
T*
CuVector<T>::data()
{
    return m_dataOnDevice;
}

template <class T>
CuVector<T>&
CuVector<T>::operator*=(const T& scalar)
{
    assertHasElements();
    OPM_CUBLAS_SAFE_CALL(detail::cublasScal(m_cuBlasHandle.get(), m_numberOfElements, &scalar, data(), 1));
    return *this;
}

template <class T>
CuVector<T>&
CuVector<T>::axpy(T alpha, const CuVector<T>& y)
{
    assertHasElements();
    assertSameSize(y);
    OPM_CUBLAS_SAFE_CALL(detail::cublasAxpy(m_cuBlasHandle.get(), m_numberOfElements, &alpha, y.data(), 1, data(), 1));
    return *this;
}

template <class T>
T
CuVector<T>::dot(const CuVector<T>& other) const
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
CuVector<T>::two_norm() const
{
    assertHasElements();
    T result = T(0);
    OPM_CUBLAS_SAFE_CALL(detail::cublasNrm2(m_cuBlasHandle.get(), m_numberOfElements, data(), 1, &result));
    return result;
}

template <typename T>
T
CuVector<T>::dot(const CuVector<T>& other, const CuVector<int>& indexSet, CuVector<T>& buffer) const
{
    return detail::innerProductAtIndices(m_dataOnDevice, other.data(), buffer.data(), indexSet.dim(), indexSet.data());
}

template <typename T>
T
CuVector<T>::two_norm(const CuVector<int>& indexSet, CuVector<T>& buffer) const
{
    // TODO: [perf] Optimize this to a single call
    return std::sqrt(this->dot(*this, indexSet, buffer));
}

template <typename T>
T
CuVector<T>::dot(const CuVector<T>& other, const CuVector<int>& indexSet) const
{
    CuVector<T> buffer(indexSet.dim());
    return detail::innerProductAtIndices(m_dataOnDevice, other.data(), buffer.data(), indexSet.dim(), indexSet.data());
}

template <typename T>
T
CuVector<T>::two_norm(const CuVector<int>& indexSet) const
{
    CuVector<T> buffer(indexSet.dim());
    // TODO: [perf] Optimize this to a single call
    return std::sqrt(this->dot(*this, indexSet, buffer));
}
template <class T>
CuVector<T>&
CuVector<T>::operator+=(const CuVector<T>& other)
{
    assertHasElements();
    assertSameSize(other);
    // TODO: [perf] Make a specialized version of this
    return axpy(1.0, other);
}

template <class T>
CuVector<T>&
CuVector<T>::operator-=(const CuVector<T>& other)
{
    assertHasElements();
    assertSameSize(other);
    // TODO: [perf] Make a specialized version of this
    return axpy(-1.0, other);
}


template <class T>
void
CuVector<T>::copyFromHost(const T* dataPointer, size_t numberOfElements)
{
    if (numberOfElements > dim()) {
        OPM_THROW(std::runtime_error,
                  fmt::format("Requesting to copy too many elements. Vector has {} elements, while {} was requested.",
                              dim(),
                              numberOfElements));
    }
    OPM_CUDA_SAFE_CALL(cudaMemcpy(data(), dataPointer, numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
void
CuVector<T>::copyToHost(T* dataPointer, size_t numberOfElements) const
{
    assertSameSize(detail::to_int(numberOfElements));
    OPM_CUDA_SAFE_CALL(cudaMemcpy(dataPointer, data(), numberOfElements * sizeof(T), cudaMemcpyDeviceToHost));
}

template <class T>
void
CuVector<T>::copyFromHost(const std::vector<T>& data)
{
    copyFromHost(data.data(), data.size());
}
template <class T>
void
CuVector<T>::copyToHost(std::vector<T>& data) const
{
    copyToHost(data.data(), data.size());
}
template class CuVector<double>;
template class CuVector<float>;
template class CuVector<int>;

} // namespace Opm::cuistl
