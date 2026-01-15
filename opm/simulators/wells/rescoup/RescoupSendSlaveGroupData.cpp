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

#include <config.h>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/simulators/wells/rescoup/RescoupSendSlaveGroupData.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>

#include <array>
#include <string>
#include <vector>
#include <tuple>

#include <fmt/format.h>

namespace Opm {

// -------------------------------------------------------
// Constructor for the RescoupTargetCalculator class
// -------------------------------------------------------
template <class Scalar, class IndexTraits>
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
RescoupSendSlaveGroupData(GroupStateHelperType& groupStateHelper)
    : groupStateHelper_{groupStateHelper}
    , reservoir_coupling_slave_{groupStateHelper.reservoirCouplingSlave()}
    , schedule_{groupStateHelper.schedule()}
    , group_state_{groupStateHelper.groupState()}
    , phase_usage_{groupStateHelper.phaseUsage()}
    , report_step_idx_{groupStateHelper.reportStepIdx()}
{
}

// ------------------
// Public methods
// ------------------

// The following data from the slave might be needed by the master process at the beginning of the
// timestep depending on the situation:
//
// 1) *Production group potentials*: If both the master group and slave groups are production groups (and
//   possibly also injection groups), and the master group has GCONPROD item 10 set to "FORM", "POTN", or "0".
//   In this case, the master group needs the slave group potentials to calculate the guide rates.
//   This could also apply to a parent of the master group, for example if the master group has a
//   fixed guide rate, and a parent group of the master group has GCONPROD item 10 set to "FORM". Then the
//   parent group needs the slave group potentials to calculate the group potential and its guide rate.
//
// 2) *Production group surface rates*: If both the master group and slave groups are production
//   groups (and possibly also injection groups), and the master group has GCONPROD item 10 set to
//   OIL, GAS, WATER, or LIQ and the phase under group control is different from the phase in GCONPROD
//   item 10. In this case, the master group needs the slave group production rates to transform guide rate
//   targets to the phase under group control. NOTE: This should also apply to parents of the master group.
//
// 3) *Production group surface rate for phase P for reinjection*: If the master group is both an injector
//   and a producer group, and the slave group is a producer group (and possibly also an injection group), and
//   the master group has GCONINJE item 3 set to REIN for a phase P (set in GCONINJE item 2). In this case,
//   the master group needs the slave group reinjection surface rates to calculate the injection target for
//   phase P. NOTE: The reinjection surface rates are essentially the same as the (surface) production rates,
//   except that if the phase P is GAS, the reinjection surface rates also includes the gas import and
//   consumption rates. NOTE: This case could also apply to a parent of the master group. NOTE: If we send
//   the surface production rates, see (2) above, then we only need to send the reinjection rate for the
//   gas phase since the other phases are not affected by the gas import and consumption rates and can be
//   computed from the surface production rates.
//
// 4) *Production group reservoir voidage replacement rates*:
//   a) If the master group is both an injector and a producer group, and the slave group is a
//     producer group (and possibly also an injection group), and the master group has
//     GCONINJE item 10 set to NETV or VOID. In this case, the master group needs the slave group
//     reservoir voidage replacement rate to calculate the guide rate. NOTE: Also applies to parents of
//     the master group.
//   b) If the master group is both an injector and a producer group, and the slave group is a
//     producer group (and possibly also an injection group), and the master group has GCONINJE item 3
//     set to VREP. In this case, the master group needs the slave group reservoir voidage replacement
//     rate to calculate the group target for the phase under group control. NOTE: Also applies to
//     parents of the master group.
//
// 5) If a parent of the master group is a pressure maintenance group (GPMAINT not supported for the
//   master group itself), it may need to know the following slave group data:
//   a) slave group reservoir voidage production rate (GPMAINT item 2 = PROD)
//   b) slave group reservoir injection rates (GPMAINT item 2 = OINJ, WINJ, or GINJ)
//   c) slave group surface injection reates (GPMAINT item 2 = OINS, WINS, or GINS)
//
// 6) If the master group or any of its parents is a sales gas control group. It will need to know
//       a) slave group surface gas production rates
//       b) slave group surface gas injection rates
//       c) slave group gas reinjection rate
//   in order to calculate the sales gas targets.
//
// 7) If a parent of the master group has RESV control mode, and the master group is available for higher level control,
//   it will need to know the slave group's reservoir production rates to calculate the master group's target rate.
//
// TODO: We could try to send only the necessary data to the master process, but that would require
//       the master process to first tell the slave group its data requirements at the current
//       report step. Which would be complicated by having to check for parent group control modes, etc.
//       Note also that even if the slave group is e.g. a pure injection group, it still might have
//       production wells, that contribute to production group at a higher level, the same goes for
//       the master group. So to be safe, we have to send all injection and production data available
//       i.e. group potentials, reservoir and surface production rates, reservoir voidage rates, reinjection
//       rate for the gas phase, and surface and reservoir injection rates.
//
template <class Scalar, class IndexTraits>
void
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
sendSlaveGroupDataToMaster()
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    if (rescoup_slave.getComm().rank() == 0) {
        this->sendSlaveGroupProductionDataToMaster_();
        this->sendSlaveGroupInjectionDataToMaster_();
    }
}


// -------------------------------------------------------
// Private methods for class RescoupSendSlaveGroupData
// -------------------------------------------------------


template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::SlaveGroupInjectionData
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupInjectionData_(std::size_t group_idx) const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
    InjectionRates surface_rates = this->createInjectionRatesFromRateVector_(
        this->group_state_.injection_surface_rates(group.name())
    );
    InjectionRates reservoir_rates = this->createInjectionRatesFromRateVector_(
        this->group_state_.injection_reservoir_rates(group.name())
    );
    return SlaveGroupInjectionData{surface_rates, reservoir_rates};
}

template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::Potentials
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupPotentials_(std::size_t group_idx) const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
    Potentials potentials;
    const auto& gr_pot = this->group_state_.get_production_group_potential(group_name);
    potentials[ReservoirCoupling::Phase::Oil] = gr_pot.oil_rate;
    potentials[ReservoirCoupling::Phase::Gas] = gr_pot.gas_rate;
    potentials[ReservoirCoupling::Phase::Water] = gr_pot.water_rate;
    return potentials;
}

template<typename Scalar, typename IndexTraits>
Scalar
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupReinjectionRateForGasPhase_(std::size_t group_idx) const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
    const int gas_phase_idx = this->phase_usage_.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
    return this->group_state_.injection_rein_rates(group.name())[gas_phase_idx];
}

template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::ProductionRates
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupSurfaceProductionRates_(std::size_t group_idx) const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
    GuideRate::RateVector production_rates = this->groupStateHelper_.getProductionGroupRateVector(group_name);
    // NOTE: GuideRate::RateVector is a vector of doubles, so we need to convert it to Scalar
    // TODO: Fix GuideRate::RateVector to be a vector of Scalars instead of doubles
    return ProductionRates{production_rates};
}

template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::ProductionRates
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupNetworkSurfaceProductionRates_(std::size_t group_idx) const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
    const auto& pu = this->phase_usage_;

    // Collect surface rates for each phase with network=true
    // This means efficiency factors are 1.0 for groups/wells with GEFAC/WEFAC item 3 = "NO"
    Scalar oil_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        oil_rate = this->groupStateHelper_.sumWellPhaseRates(
            /*res_rates=*/false, group, pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx),
            /*is_injector=*/false, /*network=*/true);
    }

    Scalar gas_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        gas_rate = this->groupStateHelper_.sumWellPhaseRates(
            /*res_rates=*/false, group, pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx),
            /*is_injector=*/false, /*network=*/true);
    }

    Scalar water_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
        water_rate = this->groupStateHelper_.sumWellPhaseRates(
            /*res_rates=*/false, group, pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx),
            /*is_injector=*/false, /*network=*/true);
    }

    // Sum across all MPI ranks since wells in a group may be owned by different ranks
    oil_rate = this->comm().sum(oil_rate);
    gas_rate = this->comm().sum(gas_rate);
    water_rate = this->comm().sum(water_rate);
    ProductionRates network_rates;
    network_rates[ReservoirCoupling::Phase::Oil] = oil_rate;
    network_rates[ReservoirCoupling::Phase::Gas] = gas_rate;
    network_rates[ReservoirCoupling::Phase::Water] = water_rate;
    return network_rates;
}

template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::ProductionRates
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupReservoirProductionRates_(std::size_t group_idx) const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
    const auto& pu = this->phase_usage_;

    // Collect reservoir rates for each phase by summing well reservoir rates
    // These rates are computed using slave's PVT properties (formation volume factors, etc.)
    // NOTE: sumWellResRates() only sums wells owned by this rank, so we need to
    // communicate the rates across all ranks using comm().sum()
    Scalar oil_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        oil_rate = this->groupStateHelper_.sumWellResRates(
            group, pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx), /*is_injector=*/false);
    }

    Scalar gas_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        gas_rate = this->groupStateHelper_.sumWellResRates(
            group, pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx), /*is_injector=*/false);
    }

    Scalar water_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
        water_rate = this->groupStateHelper_.sumWellResRates(
            group, pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx), /*is_injector=*/false);
    }

    // Sum across all MPI ranks since wells in a group may be owned by different ranks
    oil_rate = this->comm().sum(oil_rate);
    gas_rate = this->comm().sum(gas_rate);
    water_rate = this->comm().sum(water_rate);

    ProductionRates reservoir_rates;
    reservoir_rates[ReservoirCoupling::Phase::Oil] = oil_rate;
    reservoir_rates[ReservoirCoupling::Phase::Gas] = gas_rate;
    reservoir_rates[ReservoirCoupling::Phase::Water] = water_rate;
    return reservoir_rates;
}

template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::SlaveGroupProductionData
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupProductionData_(std::size_t group_idx) const
{
    SlaveGroupProductionData production_data;
    // Potentials are used to calculate guide rates for master production groups
    production_data.potentials = this->collectSlaveGroupPotentials_(group_idx);
    // Production rates are used to transform guiderate targets for master production groups
    //   from one phase to another
    production_data.surface_rates = this->collectSlaveGroupSurfaceProductionRates_(group_idx);
    // Network surface rates are used for network leaf node calculations where efficiency
    // factors may be excluded (when GEFAC/WEFAC item 3 = "NO")
    production_data.network_surface_rates =
        this->collectSlaveGroupNetworkSurfaceProductionRates_(group_idx);
    // Reservoir rates are needed when master's parent group has RESV control mode,
    // so the conversion uses slave's PVT properties rather than master's
    production_data.reservoir_rates = this->collectSlaveGroupReservoirProductionRates_(group_idx);
    production_data.voidage_rate = this->collectSlaveGroupVoidageRate_(group_idx);
    production_data.gas_reinjection_rate = this->collectSlaveGroupReinjectionRateForGasPhase_(group_idx);
    return production_data;
}

template<typename Scalar, typename IndexTraits>
Scalar
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
collectSlaveGroupVoidageRate_(std::size_t group_idx) const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    const auto& group_name = rescoup_slave.slaveGroupIdxToGroupName(group_idx);
    const Group& group = this->schedule_.getGroup(group_name, this->report_step_idx_);
    Scalar voidage_rate = this->group_state_.injection_vrep_rate(group.name());
    return voidage_rate;
}

template<typename Scalar, typename IndexTraits>
typename RescoupSendSlaveGroupData<Scalar, IndexTraits>::InjectionRates
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
createInjectionRatesFromRateVector_(const std::vector<Scalar>& rate_vector) const
{
    const auto& pu = this->phase_usage_;
    Scalar oil_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        oil_rate = rate_vector[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
    }

    Scalar gas_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        gas_rate = rate_vector[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)];
    }

    Scalar water_rate = 0.0;
    if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
        water_rate = rate_vector[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
    }
    InjectionRates injection_rates;
    // NOTE: We could here use the more concice version of the following code:
    //    return InjectionRates{{oil_rate, gas_rate, water_rate}};
    // but that would be prone to silent bugs if the order of the phases is accidentally changed.
    injection_rates[ReservoirCoupling::Phase::Oil] = oil_rate;
    injection_rates[ReservoirCoupling::Phase::Gas] = gas_rate;
    injection_rates[ReservoirCoupling::Phase::Water] = water_rate;
    return injection_rates;
}

template <class Scalar, class IndexTraits>
void
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
sendSlaveGroupProductionDataToMaster_() const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    auto num_slave_groups = rescoup_slave.numSlaveGroups();
    if (num_slave_groups == 0) {
        // History mode: no slave groups defined, nothing to send
        return;
    }
    std::vector<SlaveGroupProductionData> production_data;
    for (std::size_t group_idx = 0; group_idx < num_slave_groups; ++group_idx) {
        // NOTE: We have to send production data even if the slave group is not a production group,
        //   i.e. if it is a pure injection group or has GroupType NONE
        production_data.emplace_back(this->collectSlaveGroupProductionData_(group_idx));
    }
    rescoup_slave.sendProductionDataToMaster(production_data);
}

template <class Scalar, class IndexTraits>
void
RescoupSendSlaveGroupData<Scalar, IndexTraits>::
sendSlaveGroupInjectionDataToMaster_() const
{
    auto& rescoup_slave = this->reservoir_coupling_slave_;
    auto num_slave_groups = rescoup_slave.numSlaveGroups();
    if (num_slave_groups == 0) {
        // History mode: no slave groups defined, nothing to send
        return;
    }
    std::vector<SlaveGroupInjectionData> injection_data;
    for (std::size_t group_idx = 0; group_idx < num_slave_groups; ++group_idx) {
        // NOTE: We would like to only send injection data only if the master group is an injector,
        //   but we cannot do that because even if the master group is a pure producer, one of its parent
        //   groups could be an injector that needs the injection data.
        //if (rescoup_slave.masterGroupIsInjector(group_idx)) {//
        // NOTE: We have to send injection data even if the slave group is not an injection group,
        //   i.e. if it is a pure production group or has GroupType NONE
        injection_data.emplace_back(this->collectSlaveGroupInjectionData_(group_idx));
    }
    rescoup_slave.sendInjectionDataToMaster(injection_data);
}


template class RescoupSendSlaveGroupData<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class RescoupSendSlaveGroupData<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
