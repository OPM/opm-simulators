// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#include <config.h>
#include <opm/models/parallel/threadmanager.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

int ThreadManager::numThreads_ = 1;

void ThreadManager::registerParameters()
{
    Parameters::Register<Parameters::ThreadsPerProcess>
        ("The maximum number of threads to be instantiated per process "
         "('-1' means 'automatic')");
}

void ThreadManager::init(bool queryCommandLineParameter)
{
    if (queryCommandLineParameter) {
        numThreads_ = Parameters::Get<Parameters::ThreadsPerProcess>();

        // some safety checks. This is pretty ugly macro-magic, but so what?
#if !defined(_OPENMP)
        if (numThreads_ != 1 && numThreads_ != -1) {
            throw std::invalid_argument("OpenMP is not available. The only valid values for "
                                        "threads-per-process is 1 and -1 but it is "+std::to_string(numThreads_)+"!");
        }
        numThreads_ = 1;
#elif !defined NDEBUG && defined DUNE_INTERFACECHECK
        if (numThreads_ != 1) {
            throw std::invalid_argument("You explicitly enabled Barton-Nackman interface checking in Dune. "
                                        "The Dune implementation of this is currently incompatible with "
                                        "thread parallelism!");
        }
        numThreads_ = 1;
#else
        if (numThreads_ == 0) {
            throw std::invalid_argument("Zero threads per process are not possible: It must be at least 1, "
                                        "(or -1 for 'automatic')!");
        }
#endif

#ifdef _OPENMP
        // actually limit the number of threads
        if (numThreads_ > 0) {
            omp_set_num_threads(numThreads_);
        }
#endif
    }

#ifdef _OPENMP
    // get the number of threads which are used in the end.
    numThreads_ = omp_get_max_threads();
#endif
}

unsigned ThreadManager::threadId()
{
#ifdef _OPENMP
    return static_cast<unsigned>(omp_get_thread_num());
#else
    return 0;
#endif
}

} // namespace Opm
