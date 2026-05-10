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
}

namespace Opm {

/// \brief Reservoir-coupling flow helper bound to a BlackoilWellModel.
///
/// Holds no state of its own; all methods reach back to the parent
/// BlackoilWellModel via well_model_ and forward to the helpers in
/// opm/simulators/wells/rescoup/.  Construct once per BlackoilWellModel.
template<typename TypeTag>
class BlackoilWellModelRescoup {
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar      = GetPropType<TypeTag, Properties::Scalar>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;

public:
    explicit BlackoilWellModelRescoup(BlackoilWellModel<TypeTag>& well_model);

    /// \brief Slave-side counterpart of sendMasterGroupConstraintsToSlaves().
    ///
    /// Blocking receive of the per-master-group production targets and
    /// injection limits computed by the master.  The received values are
    /// written into the slave's group state via the receiver helper.
    /// Called from the slave's beginTimeStep first-substep handshake.
    void receiveGroupConstraintsFromMaster();

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

    /// \brief Master-side: compute per-master-group production targets and
    ///   injection limits, then dispatch them to each activated slave.
    ///
    /// Runs once per sync step.  Drives the full master-side target-
    /// distribution flow (pre-phase + Phase 1/2/3) via the constraint
    /// calculator helper.  Called from the master's beginTimeStep
    /// first-substep handshake.
    void sendMasterGroupConstraintsToSlaves();

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
    BlackoilWellModel<TypeTag>& well_model_;
};

} // namespace Opm

#include "BlackoilWellModelRescoup_impl.hpp"

#endif // RESERVOIR_COUPLING_ENABLED
#endif // OPM_BLACKOILWELLMODEL_RESCOUP_HEADER_INCLUDED
