/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_RESERVOIR_COUPLING_MASTER_REPORT_STEP_HPP
#define OPM_RESERVOIR_COUPLING_MASTER_REPORT_STEP_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

template <class Scalar>
class ReservoirCouplingMasterReportStep {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using SlaveGroupProductionData = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;

    ReservoirCouplingMasterReportStep(
        ReservoirCouplingMaster<Scalar> &master
    );
    const Parallel::Communication &comm() const { return this->master_.getComm(); }
    const std::vector<std::string>& getMasterGroupNamesForSlave(std::size_t slave_idx) const {
         return this->master_.getMasterGroupNamesForSlave(slave_idx);
    }
    const std::map<std::string, std::string>& getMasterGroupToSlaveNameMap() const {
         return this->master_.getMasterGroupToSlaveNameMap();
    }
    std::size_t getMasterGroupCanonicalIdx(
        const std::string &slave_name, const std::string &master_group_name) const;
    MPI_Comm getSlaveComm(int index) const { return this->master_.getSlaveComm(index); }
    const Potentials& getSlaveGroupPotentials(const std::string &master_group_name) const;
    void maybeReceiveGroupInfoFromSlaves();
    std::size_t numSlaveGroups(unsigned int index) const { return this->master_.numSlaveGroups(index); }
    std::size_t numSlaves() const { return this->master_.numSlavesStarted(); }
    ReservoirCoupling::Logger& logger() const { return this->master_.getLogger(); }
    void receiveInjectionDataFromSlaves();
    void receiveProductionDataFromSlaves();
    const Schedule &schedule() const { return this->master_.schedule(); }
    void sendGroupInfoToSlaves(int report_step_idx);
    void setReportStepIdx(int report_step_idx);
    bool slaveIsActivated(int index) const { return this->master_.slaveIsActivated(index); }
    const std::string &slaveName(int index) const { return this->master_.getSlaveName(index); }

private:
    void receiveGroupInfoFromSlave_(unsigned int slave_idx);
    void sendMasterGroupInjectorProducerInfoToSlave_(unsigned int slave_idx, int report_step_idx);

    ReservoirCouplingMaster<Scalar> &master_;
    // Current report step index
    int report_step_idx_;
    // Production data and potentials for each slave group
    std::map<std::string, std::vector<SlaveGroupProductionData>> slave_group_production_data_;
    // Injection data for each slave group
    std::map<std::string, std::vector<SlaveGroupInjectionData>> slave_group_injection_data_;
};
} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_MASTER_REPORT_STEP_HPP
