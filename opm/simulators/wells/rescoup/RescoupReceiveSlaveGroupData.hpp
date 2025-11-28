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

#ifndef OPM_RESCOUP_RECEIVE_SLAVE_GROUP_DATA_HPP
#define OPM_RESCOUP_RECEIVE_SLAVE_GROUP_DATA_HPP

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm {

/// @brief Receives and processes group data from slave processes in reservoir coupling
///
/// This class handles the reception of production and injection data from slave
/// processes and integrates it into the master process's group state. It is responsible for:
/// - Receiving production data (rates, potentials) from all active slave processes
/// - Receiving injection data (rates) from all active slave processes
/// - Processing and storing the data in the appropriate group state structures
/// - Coordinating with GroupStateHelper for group control calculations
///
/// The class acts as a bridge between the MPI communication layer (ReservoirCouplingMaster)
/// and the well/group management system (GroupStateHelper), ensuring that slave group
/// data is properly incorporated into the master's group control logic.
///
/// @tparam Scalar Floating-point type for rate and potential values (typically double or float)
/// @tparam IndexTraits Type traits for phase indexing
///
/// @see ReservoirCouplingMaster
/// @see GroupStateHelper
template<class Scalar, class IndexTraits>
class RescoupReceiveSlaveGroupData {
public:
    using SlaveGroupProductionData  = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using ProductionRates = ReservoirCoupling::ProductionRates<Scalar>;
    using GroupStateHelperType = GroupStateHelper<Scalar, IndexTraits>;

    /// @brief Construct a receiver for slave group data
    /// @param groupStateHelper Reference to the GroupStateHelper for accessing group state and schedule
    RescoupReceiveSlaveGroupData(GroupStateHelperType& groupStateHelper);

    /// @brief Receive and process group data from all active slave processes
    ///
    /// This method orchestrates the reception of both production and injection data
    /// from all active slave processes via the ReservoirCouplingMaster. The received
    /// data is then processed and integrated into the master's group state, making
    /// it available for group control calculations.
    ///
    /// The method performs the following operations:
    /// 1. Receives production data (rates and potentials) from all slaves
    /// 2. Receives injection data (rates) from all slaves
    /// 3. Processes and stores the data in appropriate group state structures
    ///
    /// @note This is a blocking operation that waits for all slaves to send their data
    /// @note Must be called at appropriate synchronization points in the simulation
    void receiveSlaveGroupData();
private:
    /// Reference to the GroupStateHelper for group state management
    GroupStateHelperType& groupStateHelper_;

    /// Reference to the ReservoirCouplingMaster for MPI communication with slaves
    ReservoirCouplingMaster<Scalar>& reservoir_coupling_master_;

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
#endif // OPM_RESCOUP_RECEIVE_SLAVE_GROUP_DATA_HPP
