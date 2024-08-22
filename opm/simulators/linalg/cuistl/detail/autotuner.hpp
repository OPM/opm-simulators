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
#include <opm/common/OpmLog/OpmLog.hpp>
#include <functional>
#include <utility>
#include <limits>
#include <string>

namespace Opm::cuistl::detail
{

    /// @brief Function that tests the best thread block size, assumes updating the reference will affect runtimes
    /// @tparam func function to tune
    /// @tparam ...Args types of the arguments needed to call the function
    /// @param f the function to tune, which takes the thread block size as the input
    template <typename func, typename... Args>
    int tuneThreadBlockSize(func& f, std::string descriptionOfFunction) {

        // TODO: figure out a more rigorous way of deciding how many runs will suffice?
        constexpr const int runs = 2;
        cudaEvent_t events[runs+1];

        // create the events
        for (int i = 0; i < runs + 1; ++i){
            OPM_CUDA_SAFE_CALL(cudaEventCreate(&events[i]));
        }

        // Initialize helper variables
        float bestTime = std::numeric_limits<float>::max();
        int bestBlockSize = -1;
        int interval = 64;

        // try each possible blocksize
        for (int thrBlockSize = interval; thrBlockSize <= 1024; thrBlockSize += interval){

            // record a first event, and then an event after each kernel
            OPM_CUDA_SAFE_CALL(cudaEventRecord(events[0]));
            for (int i = 0; i < runs; ++i){
                f(thrBlockSize); // runs an arbitrary function with the provided arguments
                OPM_CUDA_SAFE_CALL(cudaEventRecord(events[i+1]));
            }

            // make suret he runs are over
            OPM_CUDA_SAFE_CALL(cudaEventSynchronize(events[runs]));

            // kernel launch was valid
            if (cudaSuccess == cudaGetLastError()){
                // check if we beat the record for the fastest kernel
                for (int i = 0; i < runs; ++i){
                    float candidateBlockSizeTime;
                    OPM_CUDA_SAFE_CALL(cudaEventElapsedTime(&candidateBlockSizeTime, events[i], events[i+1]));
                    if (candidateBlockSizeTime < bestTime){
                        bestTime = candidateBlockSizeTime;
                        bestBlockSize = thrBlockSize;
                    }
                }
            }
        }

        OpmLog::info(fmt::format("{}: Tuned Blocksize: {} (fastest runtime: {}).", descriptionOfFunction, bestBlockSize, bestTime));

        return bestBlockSize;
    }

} // end namespace Opm::cuistl::detail

#endif
