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

#ifndef OPM_BLACKOILWELLMODEL_RESCOUP_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_RESCOUP_HEADER_INCLUDED

#include <opm/simulators/flow/rescoup/ReservoirCouplingEnabled.hpp>

#ifdef RESERVOIR_COUPLING_ENABLED

#include <opm/models/utils/basicproperties.hh>

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>

#include <optional>

namespace Opm {
    class DeferredLogger;
    class Schedule;
    template<class TypeTag> class BlackoilWellModel;
    template<class Scalar, class IndexTraits> class GroupStateHelper;
    template<class Scalar, class IndexTraits> class WellState;
    template<class Scalar> class ReservoirCouplingMaster;
    template<class Scalar> class ReservoirCouplingSlave;
}

namespace Opm {

/// \brief Reservoir-coupling flow helper bound to a BlackoilWellModel.
///
/// All methods reach back to the parent BlackoilWellModel via well_model_
/// (and a few cached references to commonly-accessed sub-objects) and
/// forward to the helpers in opm/simulators/wells/rescoup/.  Construct
/// once per BlackoilWellModel.
///
/// The mode-query, accessor and reservoir-coupling-facade wrappers in
/// the public section are pure forwarders -- they exist so the method
/// bodies in BlackoilWellModelRescoup_impl.hpp read as `this->X()`
/// rather than `this->well_model_.X()`.
template<typename TypeTag>
class BlackoilWellModelRescoup {
    using FluidSystem     = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar          = GetPropType<TypeTag, Properties::Scalar>;
    using IndexTraits     = typename FluidSystem::IndexTraitsType;
    using Simulator       = GetPropType<TypeTag, Properties::Simulator>;
    using ModelParameters = BlackoilModelParameters<Scalar>;
    using WellInterfacePtr = typename BlackoilWellModel<TypeTag>::WellInterfacePtr;

public:
    explicit BlackoilWellModelRescoup(BlackoilWellModel<TypeTag>& well_model);

    // === Mode queries and rescoup-facade forwarders ===

    bool isReservoirCouplingMaster() const { return well_model_.isReservoirCouplingMaster(); }
    bool isReservoirCouplingSlave()  const { return well_model_.isReservoirCouplingSlave(); }

    ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster()
    { return well_model_.reservoirCouplingMaster(); }
    const ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster() const
    { return well_model_.reservoirCouplingMaster(); }

    ReservoirCouplingSlave<Scalar>& reservoirCouplingSlave()
    { return well_model_.reservoirCouplingSlave(); }
    const ReservoirCouplingSlave<Scalar>& reservoirCouplingSlave() const
    { return well_model_.reservoirCouplingSlave(); }

    // === Common-state accessor forwarders ===

    GroupStateHelper<Scalar, IndexTraits>& groupStateHelper()
    { return well_model_.groupStateHelper(); }
    const GroupStateHelper<Scalar, IndexTraits>& groupStateHelper() const
    { return well_model_.groupStateHelper(); }

    WellState<Scalar, IndexTraits>& wellState()
    { return well_model_.wellState(); }
    const WellState<Scalar, IndexTraits>& wellState() const
    { return well_model_.wellState(); }

    const Schedule& schedule() const { return well_model_.schedule(); }

    std::vector<WellInterfacePtr>& wellContainer()
    { return well_model_.wellContainer(); }

    // === Rescoup flow methods ===

    /// \brief True if the most recent sendMasterGroupNodePressuresToSlaves
    ///   carried is_final = true (or if no send has happened yet).  Used
    ///   to gate redundant sends and to decide whether a trailing-final
    ///   send is needed when the master's outer loop exits.  Default true:
    ///   the slave is in the "not waiting for master" state until the
    ///   first non-final send happens.
    bool lastSentMasterGroupNodePressuresIsFinal() const
    { return last_sent_master_group_node_pressures_is_final_; }

    /// \brief True if at least one master group is a leaf node in the
    ///   master's extended network.  False when (a) the master is not
    ///   in rescoup-master mode, (b) the master has no active extended
    ///   network, or (c) the master's network has leaves but none of
    ///   them correspond to a master group of any slave.
    ///
    /// Used as the gate for the cross-rescoup network exchange.
    /// Queries the parsed Schedule topology so the answer is correct
    /// at every call site including the very first beginTimeStep call.
    bool masterNetworkHasMasterGroupLeaves() const;

    /// \brief Slave-side counterpart of sendCoupledNetworkActiveStatus().
    ///
    /// Blocking receive of the master's per-sync-step "is the cross-rescoup
    /// network exchange active?" flag.  Mirrored into the slave's
    /// lastReceivedMasterGroupNodePressuresIsFinal() so the slave's updateWellControlsAndNetworkIteration()
    /// gates see the master's iteration decision before updateWellControlsAndNetworkIteration() runs.
    /// Called from the slave's beginTimeStep first-substep handshake.
    void receiveCoupledNetworkActiveStatus();

    /// \brief Slave-side counterpart of sendMasterGroupConstraintsToSlaves().
    ///
    /// Blocking receive of the per-master-group production targets and
    /// injection limits computed by the master.  The received values are
    /// written into the slave's group state via the receiver helper.
    /// Called from the slave's beginTimeStep first-substep handshake.
    void receiveGroupConstraintsFromMaster();

    /// \brief Receive master-computed network-leaf node pressures and
    ///   apply them as dynamic THP limits on every slave producer
    ///   whose group has a master-supplied pressure.
    void receiveMasterGroupNodePressuresFromMaster();

    /// \brief Master-side blocking receive of per-slave production data
    ///   (rates, potentials, and injection data) from every activated slave.
    ///
    /// The data is integrated into the master's group state so it is
    /// available to the constraint calculator and summary output.
    /// Called from the master's beginTimeStep first-substep handshake and
    /// from rescoupSyncSummaryData() when a slave has fresh data to deliver.
    void receiveSlaveGroupData();

    /// \brief End-of-substep summary-data synchronisation.
    ///
    /// On the master, blocks for any pending slave production data so that
    /// the subsequent evalSummaryState() call sees up-to-date rates.  On
    /// the slave, ships production data to the master on the last substep
    /// of the sync step.  Called from timeStepSucceeded.
    void rescoupSyncSummaryData();

    /// \brief Master-side: send a single boolean to every activated slave
    ///   telling it whether the master will iterate the cross-rescoup
    ///   network exchange this sync timestep.  active = true iff at least
    ///   one master group is a leaf in the master's extended network (see
    ///   masterNetworkHasMasterGroupLeaves()).  Called from the master's
    ///   beginTimeStep first-substep handshake.
    ///
    /// Side effect: updates last_sent_master_group_node_pressures_is_final_
    /// (is_final = !active) so the master's own updateWellControlsAndNetworkIteration()
    /// gates pick up the initial state.
    void sendCoupledNetworkActiveStatus();

    /// \brief Master-side: compute per-master-group production targets and
    ///   injection limits, then dispatch them to each activated slave.
    ///
    /// Runs once per sync step.  Drives the full master-side target-
    /// distribution flow (pre-phase + Phase 1/2/3) via the constraint
    /// calculator helper.  Called from the master's beginTimeStep
    /// first-substep handshake.
    void sendMasterGroupConstraintsToSlaves();

    /// \brief Send master-computed network-leaf node pressures to each
    ///   activated slave so the slave can apply them as dynamic THP
    ///   limits on its producers in the corresponding master groups.
    /// \param is_final True if this is the master's final pressure
    ///   send for the current sync timestep; false if the master is
    ///   iterating and expects updated rates back from the slave.
    ///
    /// Side effect: updates last_sent_master_group_node_pressures_is_final_.
    void sendMasterGroupNodePressuresToSlaves(bool is_final);

    /// \brief Slave-side counterpart of receiveSlaveGroupData().
    ///
    /// Collects current production rates, potentials and injection data
    /// from the slave's local groups and ships them to the master.  Called
    /// from the slave's beginTimeStep first-substep handshake and from
    /// rescoupSyncSummaryData() on the last substep of a sync step.
    void sendSlaveGroupDataToMaster();

    /// \brief RAII guard for the deferred-logger binding.
    /// Replaces the previous BlackoilWellModel::setupRescoupScopedLogger().
    std::optional<ReservoirCoupling::ScopedLoggerGuard> setupScopedLogger(DeferredLogger& local_logger);

private:
    /// \brief Per-slave variant of masterNetworkHasMasterGroupLeaves():
    ///   true iff at least one of the given slave's master groups is a
    ///   leaf node in the master's extended network.  This is what
    ///   determines whether *that* slave participates in the cross-rescoup
    ///   network exchange; the (public) global version is the OR over all
    ///   activated slaves and gates the master's own iteration.
    bool masterNetworkHasMasterGroupLeavesForSlave_(std::size_t slave_idx) const;

    BlackoilWellModel<TypeTag>& well_model_;
    BlackoilWellModelNetwork<TypeTag>& network_;
    Simulator& simulator_;
    const ModelParameters& param_;

    // See lastSentMasterGroupNodePressuresIsFinal().
    bool last_sent_master_group_node_pressures_is_final_{true};
};

} // namespace Opm

#include "BlackoilWellModelRescoup_impl.hpp"

#endif // RESERVOIR_COUPLING_ENABLED
#endif // OPM_BLACKOILWELLMODEL_RESCOUP_HEADER_INCLUDED
