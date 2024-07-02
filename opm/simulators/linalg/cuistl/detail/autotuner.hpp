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
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <functional>
#include <utility>
#include <limits>

namespace Opm::gpuistl::detail
{

    /// @brief Function that tests the best thread block size, assumes updating the reference will affect runtimes
    /// @tparam func function to tune
    /// @tparam ...Args types of the arguments needed to call the function
    /// @param f the function to tune
    /// @param threadBlockSize reference to the thread block size that will affect kernel executions
    /// @param ...args arguments needed by the function
    template <typename func, typename... Args>
    void tuneThreadBlockSize(func f, int& threadBlockSize, Args&&... args) {

        // decide on a number of calibration runs and allocate space for the events
        constexpr const int runs = 2;
        cudaEvent_t events[runs+1];

        // create the events
        for (int i = 0; i < runs + 1; ++i){
            OPM_GPU_SAFE_CALL(cudaEventCreate(&events[i]));
        }

        // Initialize helper variables
        float bestTime = std::numeric_limits<float>::max();
        int bestBlockSize = -1;
        int interval = 64;

        // try each possible blocksize
        for (int thrBlockSize = interval; thrBlockSize <= 1024; thrBlockSize += interval){
            // update the blocksize
            threadBlockSize = thrBlockSize;

            // record a first event, and then an event after each kernel
            OPM_GPU_SAFE_CALL(cudaEventRecord(events[0]));
            for (int i = 0; i < runs; ++i){
                f(std::forward<Args>(args)...); // runs an arbitrary function with the provided arguments
                OPM_GPU_SAFE_CALL(cudaEventRecord(events[i+1]));
            }

            // make suret he runs are over
            OPM_GPU_SAFE_CALL(cudaEventSynchronize(events[runs]));

            // kernel launch was valid
            if (cudaSuccess == cudaGetLastError()){
                // check if we beat the record for the fastest kernel
                for (int i = 0; i < runs; ++i){
                    float candidateBlockSizeTime;
                    OPM_GPU_SAFE_CALL(cudaEventElapsedTime(&candidateBlockSizeTime, events[i], events[i+1]));
                    if (candidateBlockSizeTime < bestTime){
                        bestTime = candidateBlockSizeTime;
                        bestBlockSize = thrBlockSize;
                    }
                }
            }
        }
        printf("best size: %d, best time %f\n", bestBlockSize, bestTime);

        threadBlockSize = bestBlockSize;
    }

} // end namespace Opm::gpuistl::detail

#endif
