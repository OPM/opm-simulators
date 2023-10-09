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
#include <config.h>
#include <cuda_runtime.h>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/set_device.hpp>
namespace Opm::cuistl
{
void
setDevice(int mpiRank, [[maybe_unused]] int numberOfMpiRanks)
{

    int deviceCount = -1;
    cudaGetDeviceCount(&deviceCount);

    if (deviceCount <= 0) {
        // If they have CUDA enabled (ie. using a component that needs CUDA, eg. cubicgstab or CUILU0), this will fail
        // later down the line. At this point in the simulator, we can not determine if CUDA is enabled, so we can only
        // issue a warning.
        OpmLog::warning("Could not find any CUDA devices.");
        return;
    }

    // Now do a round robin kind of assignment
    // TODO: We need to be more sophistacted here. We have no guarantee this will pick the correct device.
    const auto deviceId = mpiRank % deviceCount;
    OPM_CUDA_SAFE_CALL(cudaDeviceReset());
    OPM_CUDA_SAFE_CALL(cudaSetDevice(deviceId));
    OpmLog::info("Set CUDA device to " + std::to_string(deviceId) + " (out of " + std::to_string(deviceCount)
                 + " devices).");
}
} // namespace Opm::cuistl
