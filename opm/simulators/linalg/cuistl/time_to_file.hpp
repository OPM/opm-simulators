/*
  Copyright SINTEF AS 2022

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
#ifndef OPM_CUISTL_TIME_TO_FILE_HEADER
#define OPM_CUISTL_TIME_TO_FILE_HEADER
#include <chrono>
#include <cuda_runtime.h>
#include <fstream>
#include <mpi.h>
#include <string>

#ifdef USE_INTRUSIVE_CUDA_TIMING
#define OPM_TIME_TO_FILE(basename, numberOfNonzeroes)                                                                  \
    ::Opm::cuistl::TimeToFile timer##basename(#basename, numberOfNonzeroes)
#define OPM_CU_TIME_TO_FILE(basename, numberOfNonzeroes)                                                               \
    ::Opm::cuistl::TimeToFile timer##basename(#basename, numberOfNonzeroes, true)
#else
#define OPM_TIME_TO_FILE(b, n)
#define OPM_CU_TIME_TO_FILE(b, n)
#endif

namespace
{
std::string
makeFilename(const std::string& basename)
{
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    if (worldSize == 1) {
        return basename + "_runtimes.txt";
    }

    // Number of current process
    int processId;
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    return basename + "_runtimes_" + std::to_string(processId) + ".txt";
}
} // namespace

namespace Opm::cuistl
{
struct TimeToFile {
    TimeToFile(const std::string& filename, int numberOfNonzeroes, bool doDeviceSynchronize = false)
        : filename(makeFilename(filename))
        , numberOfNonzeroes(numberOfNonzeroes)
        , start(std::chrono::high_resolution_clock::now())
        , doDeviceSynchronize(doDeviceSynchronize)
    {
    }

    ~TimeToFile()
    {
        if (doDeviceSynchronize) {
            cudaDeviceSynchronize();
        }
        const auto stop = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        std::ofstream outfile(filename, std::ios::app);

        outfile << numberOfNonzeroes << " " << duration.count() << std::endl;
    }

private:
    const std::string filename;
    const int numberOfNonzeroes;
    const std::chrono::time_point<std::chrono::high_resolution_clock> start;
    const bool doDeviceSynchronize;
};
} // namespace Opm::cuistl
#endif
