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

#include <config.h>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#if HAVE_CUDA
#include <cuda_runtime.h>
#include <cuda.h>
#include <opm/simulators/linalg/gpuistl/set_device.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#endif

namespace Opm::gpuistl {

    /*
        * Print the device name and compute capability on every rank

        If you have an AMD GPU and you have an AMD CPU you might run
        into problems with this code when using multiple MPI ranks.
        The simulation might hang because the integrated GPU in the CPU
        is detected has Radeon compute units, but it does not support ROCM.
        This is fixable my making only the GPUS on your system visible with
        ROCR_VISIBLE_DEVICES environment variable.
    */
    void printDevice()
    {
        int mpiRank = 0;
#if HAVE_CUDA
#if HAVE_MPI
        mpiRank = FlowGenericVanguard::comm().rank();
#endif

        int deviceCount = -1;
        OPM_GPU_WARN_IF_ERROR(cudaGetDeviceCount(&deviceCount));

        const auto deviceId = mpiRank % deviceCount;

        struct cudaDeviceProp props;
        OPM_GPU_WARN_IF_ERROR(cudaGetDeviceProperties(&props, deviceId));

        std::string out;
        out = fmt::format("rank: {}, GPU: {}, Compute Capability: {}.{} (device {} out of {})\n",
            mpiRank, props.name, props.major, props.minor, deviceId, deviceCount);
        auto deferred_logger = ::Opm::DeferredLogger();
        deferred_logger.info(out);

        DeferredLogger global = gatherDeferredLogger(deferred_logger, FlowGenericVanguard::comm());
        if (mpiRank == 0) {
            global.logMessages();
        }

#endif
    }

    void setDevice()
    {
#if HAVE_CUDA
#if HAVE_MPI
        Opm::gpuistl::setDevice(FlowGenericVanguard::comm().rank(), FlowGenericVanguard::comm().size());
#else
        Opm::gpuistl::setDevice(0, 1);
#endif
#endif
    }

} // namespace Opm::gpuistl
