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
#ifndef OPM_AUTOTUNER_HPP
#define OPM_AUTOTUNER_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include <functional>
#include <limits>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_resources.hpp>
#include <string>
#include <utility>

namespace Opm::gpuistl::detail
{

/// @brief Function that tests the best thread block size, assumes the provided function depends on threadblock-size
/// @tparam The type of the function to tune
/// @param f the function to tune, which takes the thread block size as the input
/// @param descriptionOfFunction Description of function
template <typename func>
int
tuneThreadBlockSize(func& f, std::string descriptionOfFunction)
{
    // This threadblock-tuner is very simple, it tests all valid block sizes divisble by 64
    // 64 is chosen so it is a multiple of the AMD wavefront size.
    // The maximum size of a threadblock is 1024, so an exhaustive search here will not be expensive
    // We time the kernel with each possible threadblock-size, and return the one
    // that gave the fastest invidivual run.

    // TODO: figure out a more rigorous way of deciding how many runs will suffice?
    constexpr const int runs = 2;
    std::array<GPUEvent, runs+1> events;

    // Initialize helper variables
    float bestTime = std::numeric_limits<float>::max();
    int bestBlockSize = -1;
    int interval = 64;

    // try each possible blocksize
    for (int thrBlockSize = interval; thrBlockSize <= 1024; thrBlockSize += interval) {

        // record a first event, and then an event after each kernel
        OPM_GPU_SAFE_CALL(cudaEventRecord(events[0].get()));
        for (int i = 0; i < runs; ++i) {
            f(thrBlockSize); // runs an arbitrary function with the provided arguments
            OPM_GPU_SAFE_CALL(cudaEventRecord(events[i + 1].get()));
        }

        // make sure the runs are over
        OPM_GPU_SAFE_CALL(cudaEventSynchronize(events[runs].get()));

        // kernel launch was valid
        if (cudaSuccess == cudaGetLastError()) {
            // check if we beat the record for the fastest kernel
            for (int i = 0; i < runs; ++i) {
                float candidateBlockSizeTime;
                OPM_GPU_SAFE_CALL(cudaEventElapsedTime(&candidateBlockSizeTime, events[i].get(), events[i + 1].get()));
                if (candidateBlockSizeTime < bestTime) { // checks if this configuration beat the current best
                    bestTime = candidateBlockSizeTime;
                    bestBlockSize = thrBlockSize;
                }
            }
        }
    }

    OpmLog::info(
        fmt::format("[Kernel tuning completed] {}: Tuned Blocksize = {}, Fastest Runtime = {}ms.", descriptionOfFunction, bestBlockSize, bestTime));

    return bestBlockSize;
}

} // end namespace Opm::gpuistl::detail

#endif
