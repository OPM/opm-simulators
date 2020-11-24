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
#include <opm/common/ErrorMacros.hpp>

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
      isOwner_(true), rankWithFirstPerf_(0),
      comm_(new Communication(Dune::MPIHelper::getLocalCommunicator())),
      commAbove_(new CommunicateAbove(*comm_))
    {}

ParallelWellInfo::ParallelWellInfo(const std::pair<std::string,bool>& well_info,
                                   Communication allComm)
    : name_(well_info.first), hasLocalCells_(well_info.second),
      rankWithFirstPerf_(-1)
{
#if HAVE_MPI
    MPI_Comm newComm;
    int color = hasLocalCells_ ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(allComm, color, allComm.rank(), &newComm);
    comm_.reset(new Communication(newComm));
#else
    comm_.reset(new Communication(Dune::MPIHelper::getLocalCommunicator()));
    rankWithFirstPerf_ = 0;
#endif
    isOwner_ = (comm_->rank() == 0);
}

void ParallelWellInfo::communicateFirstPerforation(bool hasFirst)
{
    int first = hasFirst;
    std::vector<int> firstVec(comm_->size());
    comm_->allgather(&first, 1, firstVec.data());
    auto found = std::find_if(firstVec.begin(), firstVec.end(),
                              [](int i) -> bool{ return i;});
    rankWithFirstPerf_ = found - firstVec.begin();
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

CheckDistributedWellConnections::CheckDistributedWellConnections(const Well& well,
                                                               const ParallelWellInfo& info)
    : well_(well), pwinfo_(info)
{
    foundConnections_.resize(well.getConnections().size(), 0);
}

void
CheckDistributedWellConnections::connectionFound(std::size_t i)
{
    foundConnections_[i] = 1;
}

bool
CheckDistributedWellConnections::checkAllConnectionsFound()
{
    // Ecl does not hold any information of remote connections.
    assert(pwinfo_.communication().max(foundConnections_.size()) == foundConnections_.size());
    pwinfo_.communication().sum(foundConnections_.data(),
                                foundConnections_.size());

    std::string msg = std::string("Cells with these i,j,k indices were not found ")
        + "in grid (well = " + pwinfo_.name() + "):";
    bool missingCells = false;
    auto start = foundConnections_.begin();
    for(auto conFound = start; conFound != foundConnections_.end(); ++conFound)
    {
        if (*conFound == 0)
        {
            const auto& completion = well_.getConnections()[conFound - start];
            msg = msg + " " + std::to_string(completion.getI()) + "," +
                std::to_string(completion.getJ()) + ","
                + std::to_string(completion.getK()) + " ";
            missingCells = true;
        }
    }
    if (missingCells && pwinfo_.isOwner())
    {
        OPM_THROW(std::runtime_error, msg);
    }
    return !missingCells;
}
} // end namespace Opm
