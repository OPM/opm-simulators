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

#ifndef OPM_RESERVOIR_COUPLING_SLAVE_REPORT_STEP_HPP
#define OPM_RESERVOIR_COUPLING_SLAVE_REPORT_STEP_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>

#include <map>
#include <string>
#include <utility>

namespace Opm {

// Avoid including the complete definition of ReservoirCouplingSlave here to avoid circular dependency.
template <class Scalar> class ReservoirCouplingSlave;

/// @brief Manages slave-side reservoir coupling operations for a single report step
///
/// This class encapsulates the slave process's communication with the master process
/// during a single report step in reservoir coupling simulations. It handles:
/// - Sending production data (rates, potentials) to the master process via MPI
/// - Sending injection data (rates) to the master process via MPI
/// - Coordinating data exchange between slave and master processes
///
/// The class serves as a helper to ReservoirCouplingSlave, separating the
/// report-step-specific communication logic from the overall coupling lifecycle
/// management. This separation improves code organization and makes the coupling
/// logic easier to understand and maintain.
///
/// @tparam Scalar Floating-point type for rate and potential values (typically double or float)
///
/// @note This class holds a reference to the parent ReservoirCouplingSlave object
///       and should only be used within the scope of that object's lifetime
/// @see ReservoirCouplingSlave
template <class Scalar>
class ReservoirCouplingSlaveReportStep {
public:
    using InjectionGroupTarget = ReservoirCoupling::InjectionGroupTarget<Scalar>;
    using ProductionGroupConstraints = ReservoirCoupling::ProductionGroupConstraints<Scalar>;
    using MessageTag = ReservoirCoupling::MessageTag;
    using SlaveGroupProductionData = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;

    /// @brief Construct a report step manager for the slave process
    /// @param slave Reference to the parent ReservoirCouplingSlave object
    ReservoirCouplingSlaveReportStep(
        ReservoirCouplingSlave<Scalar> &slave
    );

    /// @brief Get the MPI communicator for intra-slave communication
    /// @return Reference to the parallel communication object
    const Parallel::Communication &comm() const { return this->slave_.getComm(); }

    /// @brief Get the MPI communicator for slave-master communication
    /// @return MPI communicator handle for communication with the master process
    MPI_Comm getSlaveMasterComm() const { return this->slave_.getMasterComm(); }

    /// @brief Per-rate-type production limits received from master hierarchy.
    /// A value of -1 means no limit defined in the hierarchy for that rate type.
    struct MasterProductionLimits {
        Scalar oil_limit{-1};
        Scalar water_limit{-1};
        Scalar gas_limit{-1};
        Scalar liquid_limit{-1};
        Scalar resv_limit{-1};
    };

    /// @brief Check if a master-imposed injection target exists for a group and phase
    /// @param gname Slave group name
    /// @param phase Injection phase (e.g., Phase::WATER, Phase::GAS)
    /// @return true if the master sent an injection target for this group/phase pair
    bool hasMasterInjectionTarget(const std::string& gname, Phase phase) const;

    /// @brief Check if master-imposed per-rate-type production limits exist for a group
    /// @param gname Slave group name
    /// @return true if the master sent per-rate-type limits for this group
    bool hasMasterProductionLimits(const std::string& gname) const;

    /// @brief Check if a master-imposed production target exists for a group
    /// @param gname Slave group name
    /// @return true if the master sent a production target for this group
    bool hasMasterProductionTarget(const std::string& gname) const;

    /// @brief Check if this is the first substep within a "sync" timestep.
    /// @details This flag is used to control reservoir coupling synchronization.
    ///          Master-slave data exchange should only happen at the start of each "sync" timestep,
    ///          not on internal substeps or at retries after convergence chops.
    /// @return true if this is the first substep of a "sync" timestep, false if not
    bool isFirstSubstepOfSyncTimestep() const { return is_first_substep_of_sync_timestep_; }

    /// @brief Get the logger for reservoir coupling operations
    /// @return Reference to the logger object for this coupling session
    ReservoirCoupling::Logger& logger() const { return this->slave_.logger(); }

    /// @brief Get the master-imposed injection target and control mode for a group and phase
    /// @param gname Slave group name
    /// @param phase Injection phase
    /// @return Pair of (target rate, control mode)
    /// @throws std::out_of_range if no target exists for the given group/phase pair
    std::pair<Scalar, Group::InjectionCMode> masterInjectionTarget(const std::string& gname, Phase phase) const;

    /// @brief Get the master-imposed per-rate-type production limits for a group
    /// @param gname Slave group name
    /// @return Reference to the limits struct (-1 = no limit for that rate type)
    /// @throws std::out_of_range if no limits exist for the given group
    const MasterProductionLimits& masterProductionLimits(const std::string& gname) const;

    /// @brief Get the master-imposed production target and control mode for a group
    /// @param gname Slave group name
    /// @return Pair of (target rate, control mode)
    /// @throws std::out_of_range if no target exists for the given group
    std::pair<Scalar, Group::ProductionCMode> masterProductionTarget(const std::string& gname) const;

    /// @brief Receive the number of injection and production constraints from master
    /// @return Pair of (num_injection_targets, num_production_constraints)
    std::pair<std::size_t, std::size_t> receiveNumGroupConstraintsFromMaster() const;

    /// @brief Receive injection group targets from master and store them locally
    /// @param num_targets Number of injection targets to receive
    void receiveInjectionGroupTargetsFromMaster(std::size_t num_targets);

    /// @brief Receive production group constraints from master and store them locally
    /// @param num_targets Number of production constraints to receive
    void receiveProductionGroupConstraintsFromMaster(std::size_t num_targets);

    /// @brief Send production data to the master process
    ///
    /// This method sends production rates, potentials, and related data for all
    /// slave groups to the master process via MPI communication. The data is used
    /// by the master for group control calculations and coordination.
    ///
    /// @param production_data Vector of production data for each slave group
    ///
    /// @note This is a blocking operation that waits for the master to receive the data
    /// @note Must be called before the master attempts to receive production data
    void sendProductionDataToMaster(const std::vector<SlaveGroupProductionData> &production_data) const;

    /// @brief Send injection data to the master process
    ///
    /// This method sends injection rates and related data for all slave groups
    /// to the master process via MPI communication. The data is used by the master
    /// for group control calculations and coordination.
    ///
    /// @param injection_data Vector of injection data for each slave group
    ///
    /// @note This is a blocking operation that waits for the master to receive the data
    /// @note Must be called before the master attempts to receive injection data
    void sendInjectionDataToMaster(const std::vector<SlaveGroupInjectionData> &injection_data) const;

    /// @brief Set whether this is the first substep within a "sync" timestep.
    /// @param value true at start of sync timestep, false after first runSubStep_() call
    void setFirstSubstepOfSyncTimestep(bool value) { is_first_substep_of_sync_timestep_ = value; }

    /// @brief Get the name of this slave process
    /// @return Reference to the name string for this slave
    const std::string& slaveName() const { return this->slave_.getSlaveName(); }

    /// @brief Store a master-imposed injection target for a group and phase
    /// @param gname Slave group name
    /// @param phase Injection phase
    /// @param target Target injection rate
    /// @param cmode Injection control mode dictated by the master
    void setMasterInjectionTarget(const std::string& gname, Phase phase, Scalar target, Group::InjectionCMode cmode);

    /// @brief Store master-imposed per-rate-type production limits for a group
    /// @param gname Slave group name
    /// @param limits Per-rate-type limits (-1 = no limit for that rate type)
    void setMasterProductionLimits(const std::string& gname, const MasterProductionLimits& limits);

    /// @brief Store a master-imposed production target for a group
    /// @param gname Slave group name
    /// @param target Target production rate
    /// @param cmode Production control mode dictated by the master
    void setMasterProductionTarget(const std::string& gname, Scalar target, Group::ProductionCMode cmode);


private:
    /// @brief Generic helper method for sending data to the master process via MPI
    ///
    /// This template method encapsulates a common MPI send pattern used for both
    /// production and injection data. It handles:
    /// - Rank checking (only rank 0 sends)
    /// - MPI datatype retrieval via Dune::MPITraits
    /// - MPI_Send operation with specified message tag
    /// - Logging of the send operation
    ///
    /// @tparam DataType The type of data being sent (e.g., SlaveGroupProductionData)
    /// @param data Vector of data structures to send
    /// @param tag MPI message tag to identify the data type
    /// @param data_type_name Human-readable name for logging (e.g., "production data")
    ///
    /// @note This is a blocking operation that waits for the master to receive the data
    /// @note Only called from rank 0 of the slave's communicator
    template <class DataType>
    void sendDataToMaster_(
        const std::vector<DataType>& data,
        MessageTag tag,
        const std::string& data_type_name
    ) const;

    /// Reference to the parent ReservoirCouplingSlave object
    ReservoirCouplingSlave<Scalar> &slave_;
    // Flag to track if this is the first substep within a "sync" timestep.
    // Used to control reservoir coupling synchronization.
    bool is_first_substep_of_sync_timestep_{true};

    // Master-imposed targets and corresponding control modes, received from the master
    // process at the beginning of each sync timestep. Cleared and repopulated on every
    // receive cycle. Used by GroupStateHelper to override slave-local group targets.
    //
    // Key: slave group name. Value: (target rate, production control mode).
    std::map<std::string, std::pair<Scalar, Group::ProductionCMode>> master_production_targets_;
    // Key: (injection phase, slave group name). Value: (target rate, injection control mode).
    std::map<std::pair<Phase, std::string>, std::pair<Scalar, Group::InjectionCMode>> master_injection_targets_;
    // Per-rate-type production limits from master hierarchy. Key: slave group name.
    // A limit of -1 means no limit defined in the hierarchy for that rate type.
    std::map<std::string, MasterProductionLimits> master_production_limits_;
};
} // namespace Opm
#endif // OPM_RESERVOIR_COUPLING_SLAVE_REPORT_STEP_HPP
