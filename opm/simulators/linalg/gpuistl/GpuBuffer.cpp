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


#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/simulators/linalg/gpuistl/MiniMatrix.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include <fmt/core.h>

#include <algorithm>
#include <cstddef>

namespace Opm::gpuistl
{

// TODO: do we have to instantiate everything like this?

template class GpuBuffer<size_t>;
template class GpuBuffer<double>;
template class GpuBuffer<float>;
template class GpuBuffer<int>;
template class GpuBuffer<bool>;
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

using ThreeByThree = MiniMatrix<double, 9>;
using TwoByTwo = MiniMatrix<double, 4>;
using ThreeByThreeNeighborInfo = NeighborInfoStruct<ResidualNBInfoStruct<true, true ,true>, ThreeByThree>;
using TwoByTwoNeighborInfo = NeighborInfoStruct<ResidualNBInfoStruct<true, true ,true>, TwoByTwo>;

template class GpuBuffer<ThreeByThreeNeighborInfo>;

template GpuView<ThreeByThreeNeighborInfo> make_view(GpuBuffer<ThreeByThreeNeighborInfo>&);
template GpuView<const ThreeByThreeNeighborInfo> make_view(const GpuBuffer<ThreeByThreeNeighborInfo>&);
template GpuView<TwoByTwoNeighborInfo> make_view(GpuBuffer<TwoByTwoNeighborInfo>&);
template GpuView<const TwoByTwoNeighborInfo> make_view(const GpuBuffer<TwoByTwoNeighborInfo>&);

} // namespace Opm::gpuistl
