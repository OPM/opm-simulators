// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

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
#include <opm/simulators/flow/DamarisWriter.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <Damaris.h>
#include <fmt/format.h>

#include <limits>
#include <stdexcept>

namespace Opm::DamarisOutput {

int setPosition(const char* field, int64_t pos)
{
    int dam_err = damaris_set_position(field, &pos);
    return dam_err;
}

int setParameter(const char* field, int value)
{
    int dam_err = damaris_parameter_set(field, &value, sizeof(int));
    return dam_err;
}

int write(const char* field, const void* data)
{
    int dam_err = damaris_write(field, data);
    return dam_err;
}

int endIteration()
{
    int dam_err =  damaris_end_iteration();
    return dam_err;
}

int setupWritingPars(Parallel::Communication comm,
                     const int n_elements_local_grid,
                     std::vector<unsigned long long>& elements_rank_offsets)
{
    // one for each rank -- to be gathered from each client rank
    std::vector<unsigned long long> elements_rank_sizes(comm.size());
    // n_elements_local_grid should be the full model size
    const unsigned long long n_elements_local = n_elements_local_grid;

    // This gets the n_elements_local from all ranks and copies them to a std::vector of all the values on all ranks
    // (elements_rank_sizes[]).
    comm.allgather(&n_elements_local, 1, elements_rank_sizes.data());
    elements_rank_offsets[0] = 0ULL;
    // This scan makes the offsets to the start of each ranks grid section if each local grid data was concatenated (in
    // rank order)

    std::partial_sum(elements_rank_sizes.begin(),
                     std::prev(elements_rank_sizes.end()),
                     elements_rank_offsets.begin() + 1);

    // find the global/total size
    unsigned long long n_elements_global_max = elements_rank_offsets[comm.size() - 1];
    n_elements_global_max += elements_rank_sizes[comm.size() - 1]; // add the last ranks size to the already accumulated offset values

    if (comm.rank() == 0) {
        OpmLog::debug(fmt::format("In setupDamarisWritingPars(): n_elements_global_max = {}",
                                  n_elements_global_max));
    }

    // Set the paramater so that the Damaris servers can allocate the correct amount of memory for the variabe
    // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
    // ToDo: Do we need to check that local ranks are 0 based ?
    int dam_err = setParameter("n_elements_local", elements_rank_sizes[comm.rank()]);
    // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
    // ToDo: Do we need to check that n_elements_global_max will fit in a C int type (INT_MAX)
    if ( n_elements_global_max <= std::numeric_limits<int>::max() ) {
        setParameter("n_elements_total", n_elements_global_max);
    } else {
        if (comm.rank() == 0) {
            OpmLog::error(fmt::format("The size of the global array ({}) is"
                                      "greater than what a Damaris paramater type supports ({}).  ",
                                      n_elements_global_max, std::numeric_limits<int>::max() ));
        }
        OPM_THROW(std::runtime_error, "setupDamarisWritingPars() n_elements_global_max "
                                      "> std::numeric_limits<int>::max() " + std::to_string(dam_err));
    }

    return dam_err;
}


void
handleError(const int dam_err, Parallel::Communication comm, const std::string& message)
{
    // Find if some rank has encountered an error.
    const int isOk = (dam_err == DAMARIS_OK);
    const bool error = (comm.sum(isOk) != comm.size());

    if (error) {
        // Form error message on ranks that had error, and put it into a DeferredLogger.
        DeferredLogger logger;
        if (dam_err != DAMARIS_OK) {
            // Since the simulator will continue, this is a warning not an error
            // from the OPM Flow point of view.
            logger.warning("OPM_DAMARIS_ERROR",
                           fmt::format("Damaris error in {}, on rank {}, error string: {}",
                                       message,
                                       comm.rank(),
                                       damaris_error_string(dam_err)));
        }
        DeferredLogger global = gatherDeferredLogger(logger, comm);
        if (comm.rank() == 0) {
            global.logMessages();
        }
    }
}
}
