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
#include <sstream>
#include <cuda_runtime.h>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/set_device.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm::gpuistl
{
void
setDevice(int mpiRank, [[maybe_unused]] int numberOfMpiRanks)
{
    int deviceCount = -1;
    [[maybe_unused]] auto cuError = cudaGetDeviceCount(&deviceCount);

    if (deviceCount <= 0) {
        // If they have CUDA/HIP enabled (ie. using a component that needs CUDA, eg. gpubicgstab or CUILU0), this will fail
        // later down the line. At this point in the simulator, we can not determine if CUDA is enabled, so we can only
        // issue a warning.
        OpmLog::warning("Could not find any CUDA/HIP devices.");
        return;
    }

    // Now do a round robin kind of assignment
    // TODO: We need to be more sophistacted here. We have no guarantee this will pick the correct device.
    const auto deviceId = mpiRank % deviceCount;
    OPM_GPU_WARN_IF_ERROR(cudaDeviceReset());
    OPM_GPU_WARN_IF_ERROR(cudaSetDevice(deviceId));
}

void
printDevice (int mpiRank, int numberOfMpiRanks)
{
    int deviceCount = -1;
    OPM_GPU_WARN_IF_ERROR(cudaGetDeviceCount(&deviceCount));

    const auto deviceId = mpiRank % deviceCount;

    struct cudaDeviceProp props;
    OPM_GPU_WARN_IF_ERROR(cudaGetDeviceProperties(&props, deviceId));

    std::string out;
    if (numberOfMpiRanks > 0){
        out = fmt::format("GPU: {}, Compute Capability: {}.{} (device {} out of {})\n",
              props.name, props.major, props.minor, deviceId, deviceCount);
    }
    else{
        out = fmt::format("GPU: {}, Compute Capability: {}.{}\n",
              props.name, props.major, props.minor);
    }
    auto deferred_logger = Opm::DeferredLogger();
    deferred_logger.info(out);
    deferred_logger.logMessages();
}

} // namespace Opm::gpuistl
