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
#ifndef OPM_GPU_THREAD_UTILS_HPP
#define OPM_GPU_THREAD_UTILS_HPP
#include <cstddef>
#include <cuda.h>
#include <cuda_runtime.h>
namespace Opm::cuistl::detail
{
    constexpr inline size_t getThreads([[maybe_unused]] size_t numberOfRows)
    {
        return 1024;
    }

    inline size_t getBlocks(size_t numberOfRows)
    {
        const auto threads = getThreads(numberOfRows);
        return (numberOfRows + threads - 1) / threads;
    }

    // Kernel here is the function object of the cuda kernel
    template <class Kernel>
    inline int getCudaRecomendedThreadBlockSize(Kernel k, int suggestedThrBlockSize=-1)
    {
        if (suggestedThrBlockSize != -1){
            return suggestedThrBlockSize;
        }
        int blockSize;
        int tmpGridSize;
        std::ignore = cudaOccupancyMaxPotentialBlockSize(&tmpGridSize, &blockSize, k, 0, 0);
        return blockSize;
    }

    inline int getNumberOfBlocks(int wantedThreads, int threadBlockSize)
    {
        return (wantedThreads + threadBlockSize - 1) / threadBlockSize;
    }

} // namespace

#endif
