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

#ifndef OPM_RESCOUP_SEND_SLAVE_GROUP_DATA_HPP
#define OPM_RESCOUP_SEND_SLAVE_GROUP_DATA_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm {

/// @brief Collects and sends group data from slave process to master in reservoir coupling
///
/// This class handles the collection of production and injection data from the slave
/// process's group state and sends it to the master process via MPI. It is responsible for:
/// - Collecting production data (rates, potentials, voidage rates) for all slave groups
/// - Collecting injection data (rates) for all slave groups
/// - Converting group state data into the appropriate reservoir coupling data structures
/// - Sending the collected data to the master process via MPI communication
///
/// The class acts as a bridge between the well/group management system (GroupStateHelper)
/// and the MPI communication layer (ReservoirCouplingSlave), ensuring that slave group
/// data is properly formatted and transmitted to the master for group control calculations.
///
/// @tparam Scalar Floating-point type for rate and potential values (typically double or float)
/// @tparam IndexTraits Type traits for phase indexing
///
/// @see ReservoirCouplingSlave
/// @see GroupStateHelper
template<class Scalar, class IndexTraits>
class RescoupSendSlaveGroupData {
public:
    using SlaveGroupProductionData  = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using ProductionRates = ReservoirCoupling::ProductionRates<Scalar>;
    using InjectionRates = ReservoirCoupling::InjectionRates<Scalar>;
    using GroupStateHelperType = GroupStateHelper<Scalar, IndexTraits>;

    /// @brief Construct a sender for slave group data
    /// @param groupStateHelper Reference to the GroupStateHelper for accessing group state and schedule
    RescoupSendSlaveGroupData(GroupStateHelperType& groupStateHelper);

    /// @brief Get the communication object
    /// @return Reference to the communication object for MPI communication
    const Parallel::Communication& comm() const { return this->groupStateHelper_.comm(); }

    /// @brief Collect and send group data to the master process
    ///
    /// This method orchestrates the collection of both production and injection data
    /// for all slave groups and sends them to the master process via the
    /// ReservoirCouplingSlave. The data is collected from the current group state and
    /// formatted according to the reservoir coupling protocol.
    ///
    /// The method performs the following operations:
    /// 1. Collects production data (rates, potentials, voidage) for each slave group
    /// 2. Collects injection data (rates) for each slave group
    /// 3. Sends production data to the master via MPI
    /// 4. Sends injection data to the master via MPI
    ///
    /// @note This is a blocking operation that waits for the master to receive the data
    /// @note Must be called at appropriate synchronization points in the simulation
    void sendSlaveGroupDataToMaster();
private:
    /// @brief Collect production potentials for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return Potentials data structure containing oil, gas, and water potentials
    Potentials collectSlaveGroupPotentials_(std::size_t group_idx) const;

    /// @brief Collect complete injection data for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return SlaveGroupInjectionData structure containing all injection rates
    SlaveGroupInjectionData collectSlaveGroupInjectionData_(std::size_t group_idx) const;

    /// @brief Collect complete production data for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return SlaveGroupProductionData structure containing rates, potentials, and voidage
    SlaveGroupProductionData collectSlaveGroupProductionData_(std::size_t group_idx) const;

    /// @brief Collect the gas phase reinjection rate for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return Gas reinjection rate
    Scalar collectSlaveGroupReinjectionRateForGasPhase_(std::size_t group_idx) const;

    /// @brief Collect surface production rates for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return ProductionRates structure containing oil, gas, and water surface rates
    /// @note These rates are computed with network=false, so efficiency factors are applied
    ///       according to regular GEFAC/WEFAC settings
    ProductionRates collectSlaveGroupSurfaceProductionRates_(std::size_t group_idx) const;

    /// @brief Collect network surface production rates for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return ProductionRates structure containing oil, gas, and water surface rates for network
    /// @note These rates are computed with network=true, meaning efficiency factors are 1.0
    ///       for groups/wells with GEFAC/WEFAC item 3 = "NO"
    ProductionRates collectSlaveGroupNetworkSurfaceProductionRates_(std::size_t group_idx) const;

    /// @brief Collect reservoir production rates for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return ProductionRates structure containing oil, gas, and water reservoir rates
    /// @note These rates are needed when the master's parent group has RESV control mode,
    ///       so the conversion uses slave's PVT properties rather than master's
    ProductionRates collectSlaveGroupReservoirProductionRates_(std::size_t group_idx) const;

    /// @brief Collect the voidage rate for a specific slave group
    /// @param group_idx Index of the slave group
    /// @return Voidage rate (volume removed from reservoir)
    Scalar collectSlaveGroupVoidageRate_(std::size_t group_idx) const;

    /// @brief Convert a rate vector to InjectionRates structure
    /// @param rate_vector Vector of injection rates indexed by phase
    /// @return InjectionRates structure containing oil, gas, and water injection rates
    InjectionRates createInjectionRatesFromRateVector_(const std::vector<Scalar>& rate_vector) const;

    /// @brief Send collected production data to the master process
    ///
    /// This method sends the production data for all slave groups to the master
    /// via MPI communication through the ReservoirCouplingSlave.
    void sendSlaveGroupProductionDataToMaster_() const;

    /// @brief Send collected injection data to the master process
    ///
    /// This method sends the injection data for all slave groups to the master
    /// via MPI communication through the ReservoirCouplingSlave.
    void sendSlaveGroupInjectionDataToMaster_() const;

    /// Reference to the GroupStateHelper for group state management
    const GroupStateHelperType& groupStateHelper_;

    /// Reference to the ReservoirCouplingSlave for MPI communication with master
    ReservoirCouplingSlave<Scalar>& reservoir_coupling_slave_;

    /// Reference to the simulation schedule
    const Schedule& schedule_;

    /// Reference to the current group state
    const GroupState<Scalar>& group_state_;

    /// Reference to phase usage information for phase indexing
    const PhaseUsageInfo<IndexTraits>& phase_usage_;

    /// Current report step index
    const int report_step_idx_;
};

} // namespace Opm

#endif // OPM_RESCOUP_SEND_SLAVE_GROUP_DATA_HPP
