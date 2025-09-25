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
#include <opm/simulators/wells/WellGroupHelper.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm {

template<class Scalar, class IndexTraits>
class RescoupSendSlaveGroupData {
public:
    using SlaveGroupProductionData  = ReservoirCoupling::SlaveGroupProductionData<Scalar>;
    using SlaveGroupInjectionData = ReservoirCoupling::SlaveGroupInjectionData<Scalar>;
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
    using ProductionRates = ReservoirCoupling::ProductionRates<Scalar>;
    using InjectionRates = ReservoirCoupling::InjectionRates<Scalar>;
    using WellGroupHelperType = WellGroupHelper<Scalar, IndexTraits>;

    RescoupSendSlaveGroupData(WellGroupHelperType& wg_helper);

    void sendSlaveGroupDataToMaster();
private:

    Potentials collectSlaveGroupPotentials_(std::size_t group_idx) const;
    SlaveGroupInjectionData collectSlaveGroupInjectionData_(std::size_t group_idx) const;
    SlaveGroupProductionData collectSlaveGroupProductionData_(std::size_t group_idx) const;
    Scalar collectSlaveGroupReinjectionRateForGasPhase_(std::size_t group_idx) const;
    ProductionRates collectSlaveGroupSurfaceProductionRates_(std::size_t group_idx) const;
    Scalar collectSlaveGroupVoidageRate_(std::size_t group_idx) const;
    InjectionRates createInjectionRatesFromRateVector_(const std::vector<Scalar>& rate_vector) const;
    void sendSlaveGroupProductionDataToMaster_() const;
    void sendSlaveGroupInjectionDataToMaster_() const;

    const WellGroupHelperType& wg_helper_;
    ReservoirCouplingSlave<Scalar>& reservoir_coupling_slave_;
    const Schedule& schedule_;
    const GroupState<Scalar>& group_state_;
    const PhaseUsageInfo<IndexTraits>& phase_usage_;
    const int report_step_idx_;
};

} // namespace Opm

#endif // OPM_RESCOUP_SEND_SLAVE_GROUP_DATA_HPP
