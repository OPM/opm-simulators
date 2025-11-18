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
#include <opm/simulators/flow/rescoup/ReservoirCouplingMpiTraits.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <mpi.h>

#include <vector>

namespace Opm {

// Avoid including the complete definition of ReservoirCouplingMaster here to avoid circular dependency.
template <class Scalar> class ReservoirCouplingMaster;

/// @brief Manages master-side reservoir coupling operations for a single report step
///
/// This class encapsulates the master process's communication and coordination
/// with slave processes during a single report step in reservoir coupling simulations.
/// It handles:
/// - Receiving production and injection data from slave processes via MPI
/// - Sending group control information and constraints to slaves
/// - Managing slave group potentials and production/injection data
/// - Coordinating the exchange of information between master and slave groups
///
/// The class serves as a helper to ReservoirCouplingMaster, separating the
/// report-step-specific logic from the overall coupling lifecycle management.
/// This separation improves code organization and makes the coupling logic easier
/// to understand and maintain.
///
/// @tparam Scalar Floating-point type for rate and potential values (typically double or float)
///
/// @note This class holds references to the parent ReservoirCouplingMaster object
///       and should only be used within the scope of that object's lifetime
/// @see ReservoirCouplingMaster
template <class Scalar>
class ReservoirCouplingMasterReportStep {
public:
    using MessageTag = ReservoirCoupling::MessageTag;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using SlaveGroupProductionData = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;

    /// @brief Construct a report step manager for the master process
    /// @param master Reference to the parent ReservoirCouplingMaster object
    ReservoirCouplingMasterReportStep(
        ReservoirCouplingMaster<Scalar> &master
    );

    /// @brief Get the MPI communicator for master-slave communication
    /// @return Reference to the parallel communication object
    const Parallel::Communication &comm() const { return this->master_.getComm(); }

    /// @brief Get the names of master groups associated with a specific slave
    /// @param slave_idx Index of the slave process
    /// @return Vector of master group names that this slave is connected to
    const std::vector<std::string>& getMasterGroupNamesForSlave(std::size_t slave_idx) const {
         return this->master_.getMasterGroupNamesForSlave(slave_idx);
    }

    /// @brief Get the mapping from master group names to slave names
    /// @return Map from master group name to the corresponding slave name
    const std::map<std::string, std::string>& getMasterGroupToSlaveNameMap() const {
         return this->master_.getMasterGroupToSlaveNameMap();
    }

    /// @brief Get the canonical index for a master group
    /// @param slave_name Name of the slave process
    /// @param master_group_name Name of the master group
    /// @return Canonical index for the specified master group in the context of the slave
    std::size_t getMasterGroupCanonicalIdx(
        const std::string &slave_name, const std::string &master_group_name) const;

    /// @brief Get the MPI communicator for a specific slave process
    /// @param index Index of the slave process
    /// @return MPI communicator handle for communication with the specified slave
    MPI_Comm getSlaveComm(int index) const { return this->master_.getSlaveComm(index); }

    /// @brief Get the production potentials for a slave group
    /// @param master_group_name Name of the master group
    /// @return Reference to the potentials data for the specified group
    const Potentials& getSlaveGroupPotentials(const std::string &master_group_name) const;

    /// @brief Receive group information from slaves if needed
    ///
    /// This method conditionally receives group information from slave processes.
    /// The reception is triggered based on simulation state and timing requirements.
    void maybeReceiveGroupInfoFromSlaves();

    /// @brief Get the number of slave groups for a specific slave process
    /// @param index Index of the slave process
    /// @return Number of groups managed by the specified slave
    std::size_t numSlaveGroups(unsigned int index) const { return this->master_.numSlaveGroups(index); }

    /// @brief Get the total number of active slave processes
    /// @return Number of slaves that have been started
    std::size_t numSlaves() const { return this->master_.numSlavesStarted(); }

    /// @brief Get the logger for reservoir coupling operations
    /// @return Reference to the logger object for this coupling session
    ReservoirCoupling::Logger& logger() const { return this->master_.getLogger(); }

    /// @brief Receive injection data from all active slave processes
    ///
    /// This method receives injection rates and related data from each slave process
    /// via MPI communication. The data is organized by slave groups
    /// and stored in slave_group_injection_data_ for use in group control calculations.
    ///
    /// @note This is a blocking operation that waits for all slaves to send data
    /// @note Must be called after slaves have computed and sent their injection data
    void receiveInjectionDataFromSlaves();

    /// @brief Receive production data from all active slave processes
    ///
    /// This method receives production rates, potentials, and related data from each
    /// slave process via MPI communication. The data is organized by
    /// slave groups and stored in slave_group_production_data_ for use in group
    /// control calculations.
    ///
    /// @note This is a blocking operation that waits for all slaves to send data
    /// @note Must be called after slaves have computed and sent their production data
    void receiveProductionDataFromSlaves();

    /// @brief Get the simulation schedule
    /// @return Reference to the Schedule object containing well and group definitions
    const Schedule &schedule() const { return this->master_.schedule(); }

    /// @brief Send group control information to all slaves for a specific report step
    /// @param report_step_idx Index of the report step for which to send group info
    ///
    /// This method sends group control information (targets, constraints, well controls)
    /// from the master to each slave process. Each slave receives the information
    /// relevant to the groups it manages.
    void sendGroupInfoToSlaves(int report_step_idx);

    /// @brief Set the current report step index
    /// @param report_step_idx The report step index to set
    ///
    /// This updates the internal state to track which report step is currently
    /// being processed.
    void setReportStepIdx(int report_step_idx);

    /// @brief Check if a specific slave process has been activated
    /// @param index Index of the slave process
    /// @return true if the slave is activated, false otherwise
    bool slaveIsActivated(int index) const { return this->master_.slaveIsActivated(index); }

    /// @brief Get the name of a specific slave process
    /// @param index Index of the slave process
    /// @return Reference to the name string for the specified slave
    const std::string &slaveName(int index) const { return this->master_.getSlaveName(index); }

private:
    /// @brief Receive group information from a specific slave process
    /// @param slave_idx Index of the slave process to receive from
    void receiveGroupInfoFromSlave_(unsigned int slave_idx);

    /// @brief Send injector and producer information for master groups to a slave
    /// @param slave_idx Index of the slave process to send to
    /// @param report_step_idx Index of the report step
    void sendMasterGroupInjectorProducerInfoToSlave_(unsigned int slave_idx, int report_step_idx);

    /// Reference to the parent ReservoirCouplingMaster object
    ReservoirCouplingMaster<Scalar> &master_;

    /// Current report step index being processed
    int report_step_idx_;

    /// Production data and potentials for each slave group (map key: master group name)
    std::map<std::string, std::vector<SlaveGroupProductionData>> slave_group_production_data_;

    /// Injection data for each slave group (map key: master group name)
    std::map<std::string, std::vector<SlaveGroupInjectionData>> slave_group_injection_data_;
};
} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_MASTER_REPORT_STEP_HPP
