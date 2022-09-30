/*
  Copyright SINTEF AS 2022

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
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/impl/cublas_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/impl/cublas_wrapper.hpp>
#include <opm/simulators/linalg/cuistl/impl/cuda_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/impl/vector_operations.hpp>

#define CHECKSIZE(x)                                                                                                   \
    if (x.m_numberOfElements != m_numberOfElements) {                                                                  \
        OPM_THROW(std::invalid_argument,                                                                               \
                  "Given vector has " << x.m_numberOfElements << ", while we have " << m_numberOfElements);            \
    }
#define CHECKPOSITIVESIZE                                                                                              \
    if (m_numberOfElements <= 0) {                                                                                     \
        OPM_THROW(std::invalid_argument, "We have 0 elements");                                                        \
    }

namespace Opm::cuistl
{

template <class T>
CuVector<T>::CuVector(const std::vector<T>& data)
    : CuVector(data.data(), data.size())
{
}

template <class T>
CuVector<T>::CuVector(const int numberOfElements)
    : m_numberOfElements(numberOfElements)
    , m_cuBlasHandle(impl::CuBlasHandle::getInstance())
{
    OPM_CUDA_SAFE_CALL(cudaMalloc(&m_dataOnDevice, sizeof(T) * m_numberOfElements));
}

template <class T>
CuVector<T>::CuVector(const T* dataOnHost, const int numberOfElements)
    : CuVector(numberOfElements)
{

    OPM_CUDA_SAFE_CALL(cudaMemcpy(m_dataOnDevice, dataOnHost, m_numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
CuVector<T>&
CuVector<T>::operator=(T scalar)
{
    CHECKPOSITIVESIZE
    impl::setVectorValue(data(), m_numberOfElements, scalar);
    return *this;
}

template <class T>
CuVector<T>&
CuVector<T>::operator=(const CuVector<T>& other)
{
    CHECKPOSITIVESIZE
    CHECKSIZE(other)
    if (other.m_numberOfElements != m_numberOfElements) {
        OPM_THROW(std::invalid_argument, "Can only copy from vector of same size.");
    }
    OPM_CUDA_SAFE_CALL(
        cudaMemcpy(m_dataOnDevice, other.m_dataOnDevice, m_numberOfElements * sizeof(T), cudaMemcpyDeviceToDevice));
    return *this;
}

template <class T>
CuVector<T>::CuVector(const CuVector<T>& other)
    : CuVector(other.m_numberOfElements)
{
    CHECKPOSITIVESIZE
    CHECKSIZE(other)
    OPM_CUDA_SAFE_CALL(
        cudaMemcpy(m_dataOnDevice, other.m_dataOnDevice, m_numberOfElements * sizeof(T), cudaMemcpyDeviceToDevice));
}

template <class T>
CuVector<T>::~CuVector()
{
    OPM_CUDA_SAFE_CALL(cudaFree(m_dataOnDevice));
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
    return m_numberOfElements;
}

template <typename T>
std::vector<T>
CuVector<T>::asStdVector() const
{
    std::vector<T> temporary(m_numberOfElements);
    copyToHost(temporary);
    return temporary;
}

template <typename T>
void
CuVector<T>::setZeroAtIndexSet(const CuVector<int>& indexSet)
{
    impl::setZeroAtIndexSet(m_dataOnDevice, indexSet.dim(), indexSet.data());
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
    CHECKPOSITIVESIZE
    OPM_CUBLAS_SAFE_CALL(impl::cublasScal(m_cuBlasHandle.get(), m_numberOfElements, &scalar, data(), 1));
    return *this;
}

template <class T>
CuVector<T>&
CuVector<T>::axpy(T alpha, const CuVector<T>& y)
{
    CHECKPOSITIVESIZE
    CHECKSIZE(y)
    OPM_CUBLAS_SAFE_CALL(impl::cublasAxpy(m_cuBlasHandle.get(), m_numberOfElements, &alpha, y.data(), 1, data(), 1));
    return *this;
}

template <class T>
T
CuVector<T>::dot(const CuVector<T>& other) const
{
    CHECKPOSITIVESIZE
    CHECKSIZE(other)
    T result = T(0);
    OPM_CUBLAS_SAFE_CALL(
        impl::cublasDot(m_cuBlasHandle.get(), m_numberOfElements, data(), 1, other.data(), 1, &result));
    return result;
}
template <class T>
T
CuVector<T>::two_norm() const
{
    CHECKPOSITIVESIZE
    T result = T(0);
    OPM_CUBLAS_SAFE_CALL(impl::cublasNrm2(m_cuBlasHandle.get(), m_numberOfElements, data(), 1, &result));
    return result;
}

template <typename T>
T
CuVector<T>::dot(const CuVector<T>& other, const CuVector<int>& indexSet, CuVector<T>& buffer) const
{
    return impl::innerProductAtIndices(m_dataOnDevice, other.data(), buffer.data(), indexSet.dim(), indexSet.data());
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
    return impl::innerProductAtIndices(m_dataOnDevice, other.data(), buffer.data(), indexSet.dim(), indexSet.data());
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
    CHECKPOSITIVESIZE
    CHECKSIZE(other)
    // TODO: [perf] Make a specialized version of this
    return axpy(1.0, other);
}

template <class T>
CuVector<T>&
CuVector<T>::operator-=(const CuVector<T>& other)
{
    CHECKPOSITIVESIZE
    CHECKSIZE(other)
    // TODO: [perf] Make a specialized version of this
    return axpy(-1.0, other);
}


template <class T>
void
CuVector<T>::copyFromHost(const T* dataPointer, int numberOfElements)
{
    if (numberOfElements > dim()) {
        OPM_THROW(std::runtime_error,
                  "Requesting to copy too many elements. Vector has " << dim() << " elements, while "
                                                                      << numberOfElements << " was requested.");
    }
    OPM_CUDA_SAFE_CALL(cudaMemcpy(data(), dataPointer, numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
void
CuVector<T>::copyToHost(T* dataPointer, int numberOfElements) const
{
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
