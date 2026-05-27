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

#ifndef OPM_BLACKOILWELLMODEL_RESCOUP_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_RESCOUP_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_BLACKOILWELLMODEL_RESCOUP_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelRescoup.hpp>
#endif

#ifdef RESERVOIR_COUPLING_ENABLED

#include <opm/common/TimingMacros.hpp>

#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/rescoup/RescoupConstraintsCalculator.hpp>
#include <opm/simulators/wells/rescoup/RescoupReceiveGroupConstraints.hpp>
#include <opm/simulators/wells/rescoup/RescoupReceiveSlaveGroupData.hpp>
#include <opm/simulators/wells/rescoup/RescoupSendSlaveGroupData.hpp>

#include <cassert>

namespace Opm {

// Constructor
// -----------
template<typename TypeTag>
BlackoilWellModelRescoup<TypeTag>::
BlackoilWellModelRescoup(BlackoilWellModel<TypeTag>& well_model)
    : well_model_{well_model}
    , network_{well_model.network()}
    , simulator_{well_model.simulator()}
    , param_{well_model.param()}
{}

// Public methods alphabetically
// ------------------------------

template<typename TypeTag>
bool
BlackoilWellModelRescoup<TypeTag>::
masterNetworkHasMasterGroupLeaves() const
{
    // Query the parsed Schedule topology (`network.has_node(...)`) rather than the runtime
    // `node_pressures_` map.  The runtime map is empty on the very first beginTimeStep call
    // (network_.update has not yet run for this substep).
    if (!this->isReservoirCouplingMaster()) return false;
    const auto& rcm = this->reservoirCouplingMaster();
    const auto num_slaves = rcm.numSlaves();
    for (std::size_t s = 0; s < num_slaves; ++s) {
        if (!rcm.slaveIsActivated(s)) continue;
        if (this->masterNetworkHasMasterGroupLeavesForSlave_(s)) {
            return true;
        }
    }
    return false;
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
receiveCoupledNetworkActiveStatus()
{
    OPM_TIMEFUNCTION();
    assert(this->isReservoirCouplingSlave());
    this->reservoirCouplingSlave().receiveCoupledNetworkActiveStatusFromMaster();
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
receiveGroupConstraintsFromMaster()
{
    OPM_TIMEFUNCTION();
    RescoupReceiveGroupConstraints<Scalar, IndexTraits> constraint_receiver{
        this->well_model_.guideRateHandler(),
        this->groupStateHelper()
    };
    constraint_receiver.receiveGroupConstraintsFromMaster();
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
receiveMasterGroupNodePressuresFromMaster()
{
    OPM_TIMEFUNCTION();
    assert(this->isReservoirCouplingSlave());
    auto& rescoup_slave = this->reservoirCouplingSlave();
    const auto [num_pressures, _is_final] =
        rescoup_slave.receiveNumMasterGroupNodePressuresFromMaster();
    if (num_pressures > 0) {
        rescoup_slave.receiveMasterGroupNodePressuresFromMaster(num_pressures);
    }
    // Apply pressures as dynamic THP limits on every producer whose
    // group has a master-supplied pressure.  Wells in master groups that
    // are not network leaves are not touched (no entry in the map).
    // Mirrors the standard local-network apply pattern in
    // BlackoilWellModelNetworkGeneric::updatePressures(): when the well is
    // currently THP-controlled, also write the new THP into the WellState
    // directly, because setDynamicThpLimit() alone leaves the active
    // control's THP value stale and the subsequent well-solve would
    // converge against the old THP.
    //
    // TODO (follow-up PR): this THP-direct path is only correct when the
    // slave has no extended network between its slave group and the wells
    // (the case exercised here).  When the slave uses an extended network
    // with the slave group declared as a fixed-pressure node, the master-
    // supplied pressure should instead be installed as that node's terminal
    // pressure and the slave's network solver should propagate it down to
    // the wells.  The plumbing exists -- the solver already reads
    // Network::Node::terminal_pressure() when walking up branches -- but
    // surfacing the master-sent value into the solver needs a runtime
    // override path (e.g. a fixed-pressure override map on
    // BlackoilWellModelNetwork) so we do not mutate the parsed Schedule
    // network at runtime.  Until then, decks where the slave group is a
    // fixed-pressure node in the slave's extended network are not handled
    // correctly.
    const auto& pressures = rescoup_slave.masterGroupNodePressures();
    if (pressures.empty()) return;
    const auto& summary_state = this->well_model_.summaryState();
    auto& well_state = this->wellState();
    for (auto& well : this->wellContainer()) {
        if (!well->isProducer() || !well->wellEcl().predictionMode()) continue;
        const auto it = pressures.find(well->wellEcl().groupName());
        if (it == pressures.end()) continue;
        well->setDynamicThpLimit(it->second);
        auto& ws = well_state[well->indexOfWell()];
        if (ws.production_cmode == Well::ProducerCMode::THP) {
            ws.thp = well->getTHPConstraint(summary_state);
        }
    }
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
receiveSlaveGroupData()
{
    OPM_TIMEFUNCTION();
    assert(this->isReservoirCouplingMaster());
    RescoupReceiveSlaveGroupData<Scalar, IndexTraits> slave_group_data_receiver{
        this->groupStateHelper(),
    };
    slave_group_data_receiver.receiveSlaveGroupData();
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
rescoupSyncSummaryData()
{
    // Reservoir coupling: exchange production data between slaves and master.
    //
    // Master side: after its first substep, the master blocks here until all
    // slaves have completed the sync step and sent their production data.
    // This ensures evalSummaryState() (called next in endTimeStep) and all
    // subsequent master substeps have correct slave production rates.
    //
    // Slave side: on the last substep of the sync step, the slave sends its
    // production data to the master.  The master is already waiting at this
    // point (blocked on MPI_Recv from its first substep's timeStepSucceeded).
    if (this->isReservoirCouplingMaster()) {
        if (this->reservoirCouplingMaster().needsSlaveDataReceive()) {
            this->receiveSlaveGroupData();
            this->reservoirCouplingMaster().setNeedsSlaveDataReceive(false);
        }
    }
    if (this->isReservoirCouplingSlave()) {
        if (this->reservoirCouplingSlave().isLastSubstepOfSyncTimestep()) {
            this->sendSlaveGroupDataToMaster();
        }
    }
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
sendCoupledNetworkActiveStatus()
{
    OPM_TIMEFUNCTION();
    assert(this->isReservoirCouplingMaster());
    // Send to each activated slave a single bool: "are you connected to the
    // master's cross-rescoup network this sync timestep?", i.e. does this
    // slave have a master group that is a leaf in the master network.  This
    // is per-slave: a slave with no leaf in the master network does not
    // participate even if other slaves do.  Sent as a dedicated one-element
    // MPI message; the per-iteration node-pressure sends inside
    // updateWellControlsAndNetworkIteration() carry their own is_final flag
    // in the NumMasterGroupNodePressures header.
    auto& rescoup_master = this->reservoirCouplingMaster();
    const auto num_slaves = rescoup_master.numSlaves();
    bool any_connected = false;
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (rescoup_master.slaveIsActivated(slave_idx)) {
            const bool connected =
                this->masterNetworkHasMasterGroupLeavesForSlave_(slave_idx);
            any_connected = any_connected || connected;
            rescoup_master.sendCoupledNetworkActiveStatusToSlave(slave_idx, connected);
        }
    }
    // Mirror the *global* active state into the master's own is_final flag:
    // the master iterates the exchange iff at least one slave is connected.
    // (This must stay global -- a per-slave value would break the master's
    // own updateWellControlsAndNetworkIteration() gate.)
    this->last_sent_master_group_node_pressures_is_final_ = !any_connected;
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
sendMasterGroupConstraintsToSlaves()
{
    OPM_TIMEFUNCTION();
    // This function is called by the master process to send the group
    // constraints to the slaves.  The "will the master iterate cross-rescoup?"
    // flag is shipped separately by sendCoupledNetworkActiveStatus().
    RescoupConstraintsCalculator<Scalar, IndexTraits> constraints_calculator{
        this->well_model_.guideRateHandler(),
        this->groupStateHelper()
    };
    constraints_calculator.calculateMasterGroupConstraintsAndSendToSlaves();
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
sendMasterGroupNodePressuresToSlaves(bool is_final)
{
    OPM_TIMEFUNCTION();
    assert(this->isReservoirCouplingMaster());
    auto& rescoup_master = this->reservoirCouplingMaster();
    const auto& node_pressures = this->network_.nodePressures();
    const auto num_slaves = rescoup_master.numSlaves();
    for (std::size_t slave_idx = 0; slave_idx < num_slaves; ++slave_idx) {
        if (!rescoup_master.slaveIsActivated(slave_idx)) continue;
        std::vector<typename ReservoirCoupling::MasterGroupNodePressure<Scalar>> pressures;
        const auto& master_groups = rescoup_master.getMasterGroupNamesForSlave(slave_idx);
        for (std::size_t i = 0; i < master_groups.size(); ++i) {
            const auto it = node_pressures.find(master_groups[i]);
            if (it != node_pressures.end()) {
                pressures.push_back({i, it->second});
            }
        }
        rescoup_master.sendNumMasterGroupNodePressuresToSlave(
            slave_idx, pressures.size(), is_final);
        if (!pressures.empty()) {
            rescoup_master.sendMasterGroupNodePressuresToSlave(slave_idx, pressures);
        }
    }
    this->last_sent_master_group_node_pressures_is_final_ = is_final;
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
sendSlaveGroupDataToMaster()
{
    OPM_TIMEFUNCTION();
    assert(this->isReservoirCouplingSlave());
    RescoupSendSlaveGroupData<Scalar, IndexTraits> slave_group_data_sender{
        this->groupStateHelper()};
    slave_group_data_sender.sendSlaveGroupDataToMaster();
}

// Automatically manages the lifecycle of the DeferredLogger pointer
// in the reservoir coupling logger.  Ensures the logger is properly
// cleared when it goes out of scope, preventing dangling pointer issues:
//
// - The ScopedLoggerGuard constructor sets the logger pointer
// - When the guard goes out of scope, the destructor clears the pointer
// - Move semantics transfer ownership safely when returning from this function
//   - The moved-from guard is "nullified" and its destructor does nothing
//   - Only the final guard in the caller will clear the logger
template<typename TypeTag>
std::optional<ReservoirCoupling::ScopedLoggerGuard>
BlackoilWellModelRescoup<TypeTag>::
setupScopedLogger(DeferredLogger& local_logger)
{
    if (this->isReservoirCouplingMaster()) {
        return ReservoirCoupling::ScopedLoggerGuard{
            this->reservoirCouplingMaster().logger(),
            &local_logger
        };
    } else if (this->isReservoirCouplingSlave()) {
        return ReservoirCoupling::ScopedLoggerGuard{
            this->reservoirCouplingSlave().logger(),
            &local_logger
        };
    }
    return std::nullopt;
}

// Private methods alphabetically
// ------------------------------

template<typename TypeTag>
bool
BlackoilWellModelRescoup<TypeTag>::
masterNetworkHasMasterGroupLeavesForSlave_(std::size_t slave_idx) const
{
    // See masterNetworkHasMasterGroupLeaves() for why the parsed Schedule
    // topology is queried rather than the runtime node_pressures_ map.
    if (!this->isReservoirCouplingMaster()) return false;
    const int episodeIdx = this->simulator_.episodeIndex();
    const auto& network = this->schedule()[episodeIdx].network();
    if (!network.active()) return false;
    const auto& rcm = this->reservoirCouplingMaster();
    for (const auto& mg : rcm.getMasterGroupNamesForSlave(slave_idx)) {
        if (network.has_node(mg)) {
            return true;
        }
    }
    return false;
}


} // namespace Opm

#endif // RESERVOIR_COUPLING_ENABLED
#endif // OPM_BLACKOILWELLMODEL_RESCOUP_IMPL_HEADER_INCLUDED
