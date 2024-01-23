/*
  Copyright 2019 Equinor AS.

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
#include "ParallelRestart.hpp"

#if HAVE_MPI
#include <mpi.h>
#endif

#if HAVE_MPI
#include <opm/simulators/utils/MPISerializer.hpp>
#endif

#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/input/eclipse/EclipseState/Util/OrderedMap.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>

namespace Opm
{

RestartValue loadParallelRestart(const EclipseIO* eclIO,
                                 Action::State& actionState,
                                 SummaryState& summaryState,
                                 const std::vector<Opm::RestartKey>& solutionKeys,
                                 const std::vector<Opm::RestartKey>& extraKeys,
                                 Parallel::Communication comm)
{
#if HAVE_MPI
    RestartValue restartValues{};

    if (eclIO)
    {
        assert(comm.rank() == 0);
        restartValues = eclIO->loadRestart(actionState, summaryState, solutionKeys, extraKeys);
    }

    Parallel::MpiSerializer ser(comm);
    ser.broadcast(0, restartValues, summaryState);
    return restartValues;
#else
    (void) comm;
    return eclIO->loadRestart(actionState, summaryState, solutionKeys, extraKeys);
#endif
}

} // end namespace Opm
