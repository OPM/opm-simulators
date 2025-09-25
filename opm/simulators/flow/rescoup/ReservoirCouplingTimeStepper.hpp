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

#ifndef OPM_RESERVOIR_COUPLING_TIME_STEPPER_HPP
#define OPM_RESERVOIR_COUPLING_TIME_STEPPER_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

template <class Scalar>
class ReservoirCouplingTimeStepper {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using Seconds = ReservoirCoupling::Seconds;

    ReservoirCouplingTimeStepper(
        ReservoirCouplingMaster<Scalar> &master
    );
    const Parallel::Communication &comm() const { return this->master_.getComm(); }
    MPI_Comm getSlaveComm(int index) const { return this->master_.getSlaveComm(index); }
    std::size_t numSlaves() const { return this->master_.numSlavesStarted(); }
    ReservoirCoupling::Logger& logger() const { return this->master_.getLogger(); }
    double maybeChopSubStep(double suggested_timestep_original, double elapsed_time) const;
    void receiveNextReportDateFromSlaves();
    void resizeNextReportDates(int size) { this->slave_next_report_time_offsets_.resize(size); }
    const Schedule& schedule() const { return this->master_.schedule(); }
    bool slaveIsActivated(int index) const { return this->master_.slaveIsActivated(index); }
    const std::string &slaveName(int index) const { return this->master_.getSlaveName(index); }
    double slaveStartDate(int index) const { return this->master_.getSlaveStartDate(index); }
    double slaveActivationDate(int index) const { return this->master_.getSlaveActivationDate(index); }
    void sendNextTimeStepToSlaves(double timestep);
    void setSlaveNextReportTimeOffset(int index, double offset) {
        this->slave_next_report_time_offsets_[index] = offset;
   }

private:
    ReservoirCouplingMaster<Scalar> &master_;
    // Elapsed time from the beginning of the simulation
    std::vector<double> slave_next_report_time_offsets_;
    // Potentials for oil, gas, and water rates for each slave group
    std::map<std::string, std::vector<Potentials>> slave_group_potentials_;
};
} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_TIME_STEPPER_HPP
