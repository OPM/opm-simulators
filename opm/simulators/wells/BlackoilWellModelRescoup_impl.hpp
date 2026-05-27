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
{}

// Public methods alphabetically
// ------------------------------

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
receiveGroupConstraintsFromMaster()
{
    OPM_TIMEFUNCTION();
    RescoupReceiveGroupConstraints<Scalar, IndexTraits> constraint_receiver{
        this->well_model_.guideRateHandler(),
        this->well_model_.groupStateHelper()
    };
    constraint_receiver.receiveGroupConstraintsFromMaster();
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
receiveSlaveGroupData()
{
    OPM_TIMEFUNCTION();
    assert(this->well_model_.isReservoirCouplingMaster());
    RescoupReceiveSlaveGroupData<Scalar, IndexTraits> slave_group_data_receiver{
        this->well_model_.groupStateHelper(),
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
    if (this->well_model_.isReservoirCouplingMaster()) {
        if (this->well_model_.reservoirCouplingMaster().needsSlaveDataReceive()) {
            this->receiveSlaveGroupData();
            this->well_model_.reservoirCouplingMaster().setNeedsSlaveDataReceive(false);
        }
    }
    if (this->well_model_.isReservoirCouplingSlave()) {
        if (this->well_model_.reservoirCouplingSlave().isLastSubstepOfSyncTimestep()) {
            this->sendSlaveGroupDataToMaster();
        }
    }
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
sendSlaveGroupDataToMaster()
{
    OPM_TIMEFUNCTION();
    assert(this->well_model_.isReservoirCouplingSlave());
    RescoupSendSlaveGroupData<Scalar, IndexTraits> slave_group_data_sender{
        this->well_model_.groupStateHelper()};
    slave_group_data_sender.sendSlaveGroupDataToMaster();
}

template<typename TypeTag>
void
BlackoilWellModelRescoup<TypeTag>::
sendMasterGroupConstraintsToSlaves()
{
    OPM_TIMEFUNCTION();
    // This function is called by the master process to send the group constraints to the slaves.
    RescoupConstraintsCalculator<Scalar, IndexTraits> constraints_calculator{
        this->well_model_.guideRateHandler(),
        this->well_model_.groupStateHelper()
    };
    constraints_calculator.calculateMasterGroupConstraintsAndSendToSlaves();
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
    if (this->well_model_.isReservoirCouplingMaster()) {
        return ReservoirCoupling::ScopedLoggerGuard{
            this->well_model_.reservoirCouplingMaster().logger(),
            &local_logger
        };
    } else if (this->well_model_.isReservoirCouplingSlave()) {
        return ReservoirCoupling::ScopedLoggerGuard{
            this->well_model_.reservoirCouplingSlave().logger(),
            &local_logger
        };
    }
    return std::nullopt;
}

} // namespace Opm

#endif // RESERVOIR_COUPLING_ENABLED
#endif // OPM_BLACKOILWELLMODEL_RESCOUP_IMPL_HEADER_INCLUDED
