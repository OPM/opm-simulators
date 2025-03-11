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
#include <cuda.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <fmt/core.h>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

namespace Opm::gpuistl
{

template <class T>
GpuView<T>::GpuView(std::vector<T>& data)
    : GpuView(data.data(), data.size())
{
}

template <typename T>
std::vector<T>
GpuView<T>::asStdVector() const
{
    std::vector<T> temporary(m_numberOfElements);
    copyToHost(temporary);
    return temporary;
}

template <class T>
void
GpuView<T>::copyFromHost(const T* dataPointer, size_t numberOfElements)
{
    if (numberOfElements > size()) {
        OPM_THROW(std::runtime_error,
                  fmt::format("Requesting to copy too many elements. View has {} elements, while {} was requested.",
                              size(),
                              numberOfElements));
    }
    OPM_GPU_SAFE_CALL(cudaMemcpy(data(), dataPointer, numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
void
GpuView<T>::copyToHost(T* dataPointer, size_t numberOfElements) const
{
    assertSameSize(numberOfElements);
    OPM_GPU_SAFE_CALL(cudaMemcpy(dataPointer, data(), numberOfElements * sizeof(T), cudaMemcpyDeviceToHost));
}

template <class T>
void
GpuView<T>::copyFromHost(const std::vector<T>& data)
{
    copyFromHost(data.data(), data.size());
}
template <class T>
void
GpuView<T>::copyToHost(std::vector<T>& data) const
{
    copyToHost(data.data(), data.size());
}

template class GpuView<double>;
template class GpuView<float>;
template class GpuView<int>;
template class GpuView<std::array<double, 3>>;
template class GpuView<std::array<float, 3>>;
template class GpuView<std::array<double, 9>>;
template class GpuView<std::array<float, 9>>;

} // namespace Opm::gpuistl
