/*
  Copyright 2024 Equinor AS

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

#ifndef OPM_RESERVOIR_COUPLING_SLAVE_HPP
#define OPM_RESERVOIR_COUPLING_SLAVE_HPP

#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

class ReservoirCouplingSlave {
public:
    using MPI_Comm_Ptr = ReservoirCoupling::MPI_Comm_Ptr;
    using MessageTag = ReservoirCoupling::MessageTag;

    ReservoirCouplingSlave(
        const Parallel::Communication &comm, const Schedule &schedule, const SimulatorTimer &timer
    );
    bool activated() const { return activated_; }
    void maybeActivate(int report_step);
    void sendActivationDateToMasterProcess() const;
    void sendNextReportDateToMasterProcess() const;
    void sendSimulationStartDateToMasterProcess() const;
    void receiveMasterGroupNamesFromMasterProcess();
    double receiveNextTimeStepFromMaster();

private:
    void checkGrupSlavGroupNames_();
    double getGrupSlavActivationDate_() const;
    void saveMasterGroupNamesAsMap_(const std::vector<char>& group_names);

    const Parallel::Communication &comm_;
    const Schedule& schedule_;
    const SimulatorTimer &timer_;
    // MPI parent communicator for a slave process
    MPI_Comm_Ptr slave_master_comm_{nullptr};
    std::map<std::string, std::string> master_group_names_;
    bool activated_{false};
};

} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_SLAVE_HPP
