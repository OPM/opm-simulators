/*
  Copyright 2020 OPM-OP AS

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
#include <opm/simulators/wells/ParallelWellInfo.hpp>

namespace Opm
{

void ParallelWellInfo::DestroyComm::operator()(Communication* comm)
{
#if HAVE_MPI
    // Only delete custom communicators.
    bool del = comm
        && (*comm != Dune::MPIHelper::getLocalCommunicator())
        && (*comm != MPI_COMM_WORLD && *comm != MPI_COMM_NULL);

    if ( del )
    {
        // Not 100% nice but safe as comm is deleted anyway
        // We can only access a copy and no reference.
        MPI_Comm mpi_comm = *comm;
        MPI_Comm_free(&mpi_comm);
    }
#endif
    delete comm;
}

ParallelWellInfo::ParallelWellInfo(const std::string& name,
                                   bool hasLocalCells)
    : name_(name), hasLocalCells_ (hasLocalCells),
      isOwner_(true), comm_(new Communication(Dune::MPIHelper::getLocalCommunicator()))
    {}

ParallelWellInfo::ParallelWellInfo(const std::pair<std::string,bool>& well_info,
                                   Communication allComm)
    : name_(well_info.first), hasLocalCells_(well_info.second)
{
#if HAVE_MPI
    MPI_Comm newComm;
    int color = hasLocalCells_ ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(allComm, color, allComm.rank(), &newComm);
    comm_.reset(new Communication(newComm));
#else
    comm_.reset(new Communication(Dune::MPIHelper::getLocalCommunicator()));
#endif
    isOwner_ = (comm_->rank() == 0);
}

bool operator<(const ParallelWellInfo& well1, const ParallelWellInfo& well2)
{
    return well1.name() < well2.name() || (! (well2.name() < well1.name()) && well1.hasLocalCells() < well2.hasLocalCells());
}

bool operator==(const ParallelWellInfo& well1, const ParallelWellInfo& well2)
{
    bool ret = well1.name() == well2.name() && well1.hasLocalCells() == well2.hasLocalCells()
        && well1.isOwner() == well2.isOwner();
#if HAVE_MPI
    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
    ret = ret &&
        static_cast<MPIComm>(well1.communication()) == static_cast<MPIComm>(well2.communication());
#endif
    return ret;
}

bool operator!=(const ParallelWellInfo& well1, const ParallelWellInfo& well2)
{
    return ! (well1 == well2);
}

bool operator<(const std::pair<std::string, bool>& pair, const ParallelWellInfo& well)
{
    return pair.first < well.name() || ( !( well.name() < pair.first ) && pair.second < well.hasLocalCells() );
}

bool operator<( const ParallelWellInfo& well, const std::pair<std::string, bool>& pair)
{
    return well.name() < pair.first || ( !( pair.first < well.name() ) && well.hasLocalCells() < pair.second );
}

bool operator==(const std::pair<std::string, bool>& pair, const ParallelWellInfo& well)
{
    return pair.first == well.name() && pair.second == well.hasLocalCells();
}

bool operator==(const ParallelWellInfo& well, const std::pair<std::string, bool>& pair)
{
    return pair == well;
}

bool operator!=(const std::pair<std::string, bool>& pair, const ParallelWellInfo& well)
{
    return pair.first != well.name() || pair.second != well.hasLocalCells();
}

bool operator!=(const ParallelWellInfo& well, const std::pair<std::string, bool>& pair)
{
    return pair != well;
}
} // end namespace Opm
