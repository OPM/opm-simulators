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
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

// Avoid including the complete definition of ReservoirCouplingMaster here to avoid circular dependency.
template <class Scalar> class ReservoirCouplingMaster;

/// @brief Manages time stepping coordination between master and slave processes
///
/// This class handles the synchronization of time steps and report dates between
/// the master and slave processes in reservoir coupling simulations. It is responsible for:
/// - Coordinating time step sizes across master and slave processes
/// - Receiving next report dates from slave processes
/// - Sending time step information to slaves
/// - Adjusting time steps to align with slave report boundaries
/// - Managing slave activation dates and simulation start dates
///
/// The class ensures that the master and all slaves advance through simulation time
/// in a coordinated manner, respecting individual slave report schedules while
/// maintaining overall simulation consistency.
///
/// @tparam Scalar Floating-point type for time and rate values (typically double or float)
///
/// @note This class holds a reference to the parent ReservoirCouplingMaster object
///       and should only be used within the scope of that object's lifetime
/// @see ReservoirCouplingMaster
template <class Scalar>
class ReservoirCouplingTimeStepper {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using Seconds = ReservoirCoupling::Seconds;

    /// @brief Construct a time stepper for coordinating master-slave time stepping
    /// @param master Reference to the parent ReservoirCouplingMaster object
    ReservoirCouplingTimeStepper(
        ReservoirCouplingMaster<Scalar> &master
    );

    /// @brief Get the MPI communicator for master-to-master communication
    /// @return Reference to the parallel communication object
    const Parallel::Communication &comm() const { return this->master_.getComm(); }

    /// @brief Get the MPI communicator for a specific slave process
    /// @param index Index of the slave process
    /// @return MPI communicator handle for communication with the specified slave
    MPI_Comm getSlaveComm(int index) const { return this->master_.getSlaveComm(index); }

    /// @brief Get the total number of active slave processes
    /// @return Number of slaves that have been started
    std::size_t numSlaves() const { return this->master_.numSlavesStarted(); }

    /// @brief Get the logger for reservoir coupling operations
    /// @return Reference to the logger object for this coupling session
    ReservoirCoupling::Logger& logger() const { return this->master_.logger(); }

    /// @brief Potentially adjust time step to align with slave report boundaries
    ///
    /// This method checks if the suggested time step would cause the master simulation
    /// to advance past any slave's next report time. If so, it "chops" the time step
    /// to ensure the master stops at the slave's report boundary, allowing for proper
    /// synchronization between master and slave processes.
    ///
    /// @param suggested_timestep_original The original suggested time step size (in seconds)
    /// @param elapsed_time Elapsed time from the beginning of the simulation (in seconds)
    /// @return Adjusted time step that respects slave report boundaries
    ///
    /// @note This is crucial for maintaining time step synchronization in coupled simulations
    double maybeChopSubStep(double suggested_timestep_original, double elapsed_time) const;

    /// @brief Receive next report dates from all active slave processes
    ///
    /// This method receives the next report time from each slave process via MPI
    /// communication. The times are stored as offsets from the simulation start time
    /// and used to coordinate time stepping between master and slaves.
    ///
    /// @note This is a blocking operation that waits for all slaves to send their data
    /// @note Must be called after slaves have determined their next report times
    void receiveNextReportDateFromSlaves();

    /// @brief Resize the internal storage for slave next report times
    /// @param size Number of slave processes to allocate storage for
    void resizeNextReportDates(int size) { this->slave_next_report_time_offsets_.resize(size); }

    /// @brief Get the simulation schedule
    /// @return Reference to the Schedule object containing timing and control information
    const Schedule& schedule() const { return this->master_.schedule(); }

    /// @brief Check if a specific slave process has been activated
    /// @param index Index of the slave process
    /// @return true if the slave is activated, false otherwise
    bool slaveIsActivated(int index) const { return this->master_.slaveIsActivated(index); }

    /// @brief Get the name of a specific slave process
    /// @param index Index of the slave process
    /// @return Reference to the name string for the specified slave
    const std::string &slaveName(int index) const { return this->master_.getSlaveName(index); }

    /// @brief Get the simulation start date for a specific slave
    /// @param index Index of the slave process
    /// @return Simulation start date as time offset (in seconds)
    double slaveStartDate(int index) const { return this->master_.getSlaveStartDate(index); }

    /// @brief Get the activation date for a specific slave
    /// @param index Index of the slave process
    /// @return Activation date as time offset (in seconds)
    double slaveActivationDate(int index) const { return this->master_.getSlaveActivationDate(index); }

    /// @brief Send the next time step size to all active slave processes
    ///
    /// This method broadcasts the time step size that slaves should use for their
    /// next simulation step via MPI communication. This ensures all processes advance
    /// through simulation time in a coordinated manner.
    ///
    /// @param timestep Time step size to send to slaves (in seconds)
    ///
    /// @note This is a blocking operation that waits for all slaves to receive the data
    void sendNextTimeStepToSlaves(double timestep);

    /// @brief Set the next report time offset for a specific slave
    /// @param index Index of the slave process
    /// @param offset Time offset from simulation start to the slave's next report time (in seconds)
    void setSlaveNextReportTimeOffset(int index, double offset) {
        this->slave_next_report_time_offsets_[index] = offset;
   }

private:
    /// Reference to the parent ReservoirCouplingMaster object
    ReservoirCouplingMaster<Scalar> &master_;

    /// Time offsets from simulation start to each slave's next report time (in seconds)
    std::vector<double> slave_next_report_time_offsets_;

    /// Production potentials for oil, gas, and water rates for each slave group
    /// (map key: master group name)
    std::map<std::string, std::vector<Potentials>> slave_group_potentials_;
};
} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_TIME_STEPPER_HPP
