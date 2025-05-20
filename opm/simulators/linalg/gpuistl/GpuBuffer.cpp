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
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include <fmt/core.h>

#include <algorithm>
#include <cstddef>

namespace Opm::gpuistl
{

template <class T>
GpuBuffer<T>::GpuBuffer(const std::vector<T>& data)
    : GpuBuffer(data.data(), data.size())
{
}

template <class T>
GpuBuffer<T>::GpuBuffer(const size_t numberOfElements)
    : m_numberOfElements(numberOfElements)
{
    OPM_GPU_SAFE_CALL(cudaMalloc(&m_dataOnDevice, sizeof(T) * m_numberOfElements));
}

template <class T>
GpuBuffer<T>::GpuBuffer(const T* dataOnHost, const size_t numberOfElements)
    : GpuBuffer(numberOfElements)
{
    OPM_GPU_SAFE_CALL(cudaMemcpy(
        m_dataOnDevice, dataOnHost, m_numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
GpuBuffer<T>::GpuBuffer(const GpuBuffer<T>& other)
    : GpuBuffer(other.m_numberOfElements)
{
    assertSameSize(other);
    if (m_numberOfElements == 0) {
        return;
    }
    OPM_GPU_SAFE_CALL(cudaMemcpy(m_dataOnDevice,
                                  other.m_dataOnDevice,
                                  m_numberOfElements * sizeof(T),
                                  cudaMemcpyDeviceToDevice));
}

template <class T>
GpuBuffer<T>::~GpuBuffer()
{
    OPM_GPU_WARN_IF_ERROR(cudaFree(m_dataOnDevice));
}

template <typename T>
typename GpuBuffer<T>::size_type
GpuBuffer<T>::size() const
{
    return m_numberOfElements;
}

template <typename T>
void
GpuBuffer<T>::resize(size_t newSize)
{
    if (newSize < 1) {
        OPM_THROW(std::invalid_argument, "Setting a GpuBuffer size to a non-positive number is not allowed");
    }

    if (m_numberOfElements == 0) {
        // We have no data, so we can just allocate new memory
        OPM_GPU_SAFE_CALL(cudaMalloc(&m_dataOnDevice, sizeof(T) * newSize));
    }
    else {
        // Allocate memory for temporary buffer
        T* tmpBuffer = nullptr;
        OPM_GPU_SAFE_CALL(cudaMalloc(&tmpBuffer, sizeof(T) * m_numberOfElements));

        // Move the data from the old to the new buffer with truncation
        size_t sizeOfMove = std::min({m_numberOfElements, newSize});
        OPM_GPU_SAFE_CALL(cudaMemcpy(tmpBuffer,
                                    m_dataOnDevice,
                                    sizeOfMove * sizeof(T),
                                    cudaMemcpyDeviceToDevice));

        // free the old buffer
        OPM_GPU_SAFE_CALL(cudaFree(m_dataOnDevice));

        // swap the buffers
        m_dataOnDevice = tmpBuffer;
    }

    // update size
    m_numberOfElements = newSize;
}

template <typename T>
std::vector<T>
GpuBuffer<T>::asStdVector() const
{
    std::vector<T> temporary(m_numberOfElements);
    copyToHost(temporary);
    return temporary;
}

template <typename T>
void
GpuBuffer<T>::assertSameSize(const GpuBuffer<T>& x) const
{
    assertSameSize(x.m_numberOfElements);
}

template <typename T>
void
GpuBuffer<T>::assertSameSize(size_t size) const
{
    if (size != m_numberOfElements) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Given buffer has {}, while we have {}.", size, m_numberOfElements));
    }
}

template <typename T>
void
GpuBuffer<T>::assertHasElements() const
{
    if (m_numberOfElements <= 0) {
        OPM_THROW(std::invalid_argument, "We have 0 elements");
    }
}

template <typename T>
T*
GpuBuffer<T>::data()
{
    return m_dataOnDevice;
}

template <typename T>
const T*
GpuBuffer<T>::data() const
{
    return m_dataOnDevice;
}

template <class T>
void
GpuBuffer<T>::copyFromHost(const T* dataPointer, size_t numberOfElements)
{
    if (numberOfElements > size()) {
        OPM_THROW(std::runtime_error,
                  fmt::format("Requesting to copy too many elements. buffer has {} elements, while {} was requested.",
                              size(),
                              numberOfElements));
    }
    OPM_GPU_SAFE_CALL(cudaMemcpy(data(), dataPointer, numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
void
GpuBuffer<T>::copyToHost(T* dataPointer, size_t numberOfElements) const
{
    assertSameSize(numberOfElements);
    OPM_GPU_SAFE_CALL(cudaMemcpy(dataPointer, data(), numberOfElements * sizeof(T), cudaMemcpyDeviceToHost));
}

template <class T>
void
GpuBuffer<T>::copyFromHost(const std::vector<T>& data)
{
    copyFromHost(data.data(), data.size());
}

template <class T>
void
GpuBuffer<T>::copyToHost(std::vector<T>& data) const
{
    copyToHost(data.data(), data.size());
}

template class GpuBuffer<size_t>;
template class GpuBuffer<double>;
template class GpuBuffer<float>;
template class GpuBuffer<int>;
template class GpuBuffer<std::byte>;
template class GpuBuffer<std::array<double, 3>>;
template class GpuBuffer<std::array<float, 3>>;
template class GpuBuffer<std::array<double, 9>>;
template class GpuBuffer<std::array<float, 9>>;

template <class T>
GpuView<T> make_view(GpuBuffer<T>& buf) {
    return GpuView<T>(buf.data(), buf.size());
}

template GpuView<double> make_view(GpuBuffer<double>&);
template GpuView<float> make_view(GpuBuffer<float>&);
template GpuView<int> make_view(GpuBuffer<int>&);
template GpuView<std::array<double, 3>> make_view(GpuBuffer<std::array<double, 3>>&);
template GpuView<std::array<float, 3>> make_view(GpuBuffer<std::array<float, 3>>&);
template GpuView<std::array<double, 9>> make_view(GpuBuffer<std::array<double, 9>>&);
template GpuView<std::array<float, 9>> make_view(GpuBuffer<std::array<float, 9>>&);

template <class T>
GpuView<const T> make_view(const GpuBuffer<T>& buf) {
    return GpuView<const T>(buf.data(), buf.size());
}

template GpuView<const double> make_view(const GpuBuffer<double>&);
template GpuView<const float> make_view(const GpuBuffer<float>&);
template GpuView<const int> make_view(const GpuBuffer<int>&);
template GpuView<const std::array<double, 3>> make_view(const GpuBuffer<std::array<double, 3>>&);
template GpuView<const std::array<float, 3>> make_view(const GpuBuffer<std::array<float, 3>>&);
template GpuView<const std::array<double, 9>> make_view(const GpuBuffer<std::array<double, 9>>&);
template GpuView<const std::array<float, 9>> make_view(const GpuBuffer<std::array<float, 9>>&);

} // namespace Opm::gpuistl
