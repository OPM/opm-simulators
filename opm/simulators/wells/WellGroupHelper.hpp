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
#ifndef OPM_WELLGROUPHELPER_HEADER_INCLUDED
#define OPM_WELLGROUPHELPER_HEADER_INCLUDED

#include <opm/common/TimingMacros.hpp>
#include <opm/input/eclipse/Schedule/Group/GPMaint.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm {

template<typename Scalar, typename IndexTraits>
class WellGroupHelper {
public:
    // RAII guard for temporarily setting wellstate pointer
    class WellStateGuard {
    public:
        WellStateGuard(WellGroupHelper& wgHelper, WellState<Scalar, IndexTraits>& well_state)
            : wgHelper_{wgHelper}
            , previous_state_ptr_{wgHelper_.well_state_}
        {
            // Set the new state directly
            wgHelper_.well_state_ = &well_state;
        }

        ~WellStateGuard() {
            // Restore the previous state
            wgHelper_.well_state_ = previous_state_ptr_;
        }

        // Delete copy and move operations
        WellStateGuard(const WellStateGuard&) = delete;
        WellStateGuard& operator=(const WellStateGuard&) = delete;
        WellStateGuard(WellStateGuard&&) = delete;
        WellStateGuard& operator=(WellStateGuard&&) = delete;

    private:
        WellGroupHelper& wgHelper_;
        WellState<Scalar, IndexTraits>* previous_state_ptr_;
    };

    // RAII guard for temporarily setting groupstate pointer
    class GroupStateGuard {
        public:
            GroupStateGuard(WellGroupHelper& wgHelper, GroupState<Scalar>& group_state)
                : wgHelper_{wgHelper}
                , previous_state_ptr_{wgHelper_.group_state_}
            {
                // Set the new state directly
                wgHelper_.group_state_ = &group_state;
            }

            ~GroupStateGuard() {
                // Restore the previous state
                wgHelper_.group_state_ = previous_state_ptr_;
            }

            // Delete copy and move operations
            GroupStateGuard(const GroupStateGuard&) = delete;
            GroupStateGuard& operator=(const GroupStateGuard&) = delete;
            GroupStateGuard(GroupStateGuard&&) = delete;
            GroupStateGuard& operator=(GroupStateGuard&&) = delete;

        private:
            WellGroupHelper& wgHelper_;
            GroupState<Scalar>* previous_state_ptr_;
        };

    WellGroupHelper(
        const Schedule& schedule,
        const SummaryState& summary_state,
        const GuideRate& guide_rate,
        const PhaseUsageInfo<IndexTraits>& phase_usage_info
    );

    void accumulateGroupEfficiencyFactor(const Group& group, Scalar& factor);
    std::pair<bool, Scalar> checkGroupConstraintsInj(
        const std::string& name,
        const std::string& parent,
        const Group& group,
        const Scalar* rates,
        const Phase injection_phase,
        const Scalar efficiency_factor,
        const std::vector<Scalar>& resv_coeff,
        const bool check_guide_rate
    ) const;
    std::pair<bool, Scalar> checkGroupConstraintsProd(
        const std::string& name,
        const std::string& parent,
        const Group& group,
        const Scalar* rates,
        const Scalar efficiency_factor,
        const std::vector<Scalar>& resv_coeff,
        const bool check_guide_rate
    ) const;
    std::map<std::string, Scalar> computeNetworkPressures(
        const Network::ExtNetwork& network,
        const VFPProdProperties<Scalar>& vfp_prod_props,
        const Parallel::Communication& comm
    ) const;
    DeferredLogger& deferredLogger() const { return *this->deferred_logger_; }
    Scalar getGuideRate(const std::string& name, const GuideRateModel::Target target) const;
    GuideRate::RateVector getProductionGroupRateVector(const std::string& group_name) const;
    Scalar getWellGroupTargetInjector(
        const std::string& name,
        const std::string& parent,
        const Group& group,
        const Scalar* rates,
        const Phase injection_phase,
        const Scalar efficiency_factor,
        const std::vector<Scalar>& resv_coeff
    );
    Scalar getWellGroupTargetProducer(
        const std::string& name,
        const std::string& parent,
        const Group& group,
        const Scalar* rates,
        const Scalar efficiency_factor,
        const std::vector<Scalar>& resv_coeff
    );
    GuideRate::RateVector getWellRateVector(const std::string& name) const;
    std::vector<std::string> groupChainTopBot(const std::string& bottom, const std::string& top) const;
    int groupControlledWells(
        const std::string& group_name,
        const std::string& always_included_child,
        const bool is_production_group,
        const Phase injection_phase
    ) const;
    const GroupState<Scalar>& groupState() const { return *this->group_state_; }
    GroupState<Scalar>& groupState() { return *this->group_state_; }
    const PhaseUsageInfo<IndexTraits>& phaseUsageInfo() const { return this->phase_usage_info_; }
    GroupStateGuard pushGroupState(GroupState<Scalar>& group_state) {
        return GroupStateGuard(*this, group_state);
    }
    WellStateGuard pushWellState(WellState<Scalar, IndexTraits>& well_state) {
        return WellStateGuard(*this, well_state);
    }
    void setCmodeGroup(const Group& group);
    void setLogger(DeferredLogger* deferred_logger) { deferred_logger_ = deferred_logger; }
    template <class AverageRegionalPressureType>
    void setRegionAveragePressureCalculator(
        const Group& group,
        const FieldPropsManager& fp,
        std::map<std::string, std::unique_ptr<AverageRegionalPressureType>>& regional_average_pressure_calculator
    );
    Scalar sumSolventRates(const Group& group, const bool is_injector);
    Scalar sumWellResRates(const Group& group, const int phase_pos, const bool injector) const;
    Scalar sumWellSurfaceRates(const Group& group, const int phase_pos, const bool injector) const;
    Scalar sumWellPhaseRates(
        bool res_rates,
        const Group& group,
        const int phase_pos,
        const bool injector,
        const bool network = false
    ) const;
    template <class RegionalValues> void updateGpMaintTargetForGroups(
        const Group& group, const RegionalValues& regional_values, const double dt
    );
    int updateGroupControlledWells(const bool is_production_group, const Phase injection_phase);
    void updateGroupProductionRates(const Group& group);
    void updateGroupTargetReduction(const Group& group, const bool is_injector);
    void updateNetworkLeafNodeProductionRates();
    void updateREINForGroups(const Group& group, const bool sum_rank);
    void updateReservoirRatesInjectionGroups(const Group& group);
    void updateVREPForGroups(const Group& group);
    void updateState(
        WellState<Scalar, IndexTraits>& well_state, GroupState<Scalar>& group_state, int report_step
    );
    void updateSurfaceRatesInjectionGroups(const Group& group);
    void updateWellRates(const Group& group, const WellState<Scalar, IndexTraits>& well_state_nupcol);
    const WellState<Scalar, IndexTraits>& wellState() const { return *this->well_state_; }
    WellState<Scalar, IndexTraits>& wellState() { return *this->well_state_; }
    void updateWellRatesFromGroupTargetScale(const Scalar scale, const Group& group, bool is_injector);
    std::pair<std::optional<std::string>, Scalar> worstOffendingWell(
        const Group& group,
        const Group::ProductionCMode& offended_control,
        const Parallel::Communication& comm
    ) const;
private:
    std::string controlGroup_(const Group& group);
    GuideRate::RateVector getGuideRateVector_(const std::vector<Scalar>& rates) const;
    bool isInGroupChainTopBot_(const std::string& bottom, const std::string& top) const;
    int phaseToActivePhaseIdx_(const Phase phase);
    Scalar satelliteInjectionRate_(
        const ScheduleState& sched,
        const Group& group,
        const int phase_pos,
        bool res_rates
    ) const;
    Scalar satelliteProductionRate_(
        const ScheduleState& sched,
        const Group& group,
        const GSatProd::GSatProdGroupProp::Rate rate_comp,
        bool res_rates
    ) const;
    std::optional<GSatProd::GSatProdGroupProp::Rate> selectRateComponent_(const int phase_pos) const;
    int updateGroupControlledWellsRecursive_(
        const std::string& group_name, const bool is_production_group, const Phase injection_phase);
    void updateGroupTargetReductionRecursive_(
        const Group& group, const bool is_injector, std::vector<Scalar>& group_target_reduction);

    const Schedule& schedule_;
    const SummaryState& summary_state_;
    WellState<Scalar, IndexTraits>* well_state_{nullptr};
    GroupState<Scalar>* group_state_{nullptr};
    const GuideRate& guide_rate_;
    int report_step_{0};
    DeferredLogger* deferred_logger_{nullptr};
    // NOTE: The phase usage info seems to be read-only throughout the simulation, so it should be safe
    // to store a reference to it here.
    const PhaseUsageInfo<IndexTraits>& phase_usage_info_;
};

// -----------------------------------------------------------------------------
// Template member function implementations
// -----------------------------------------------------------------------------

// NOTE: This template member function is defined here in the header because the
//       AverageRegionalPressureType type parameter depends on derived class types that are
//       not available when WellGroupHelper.cpp is compiled. Template functions
//       must be visible at their instantiation point.
//       See WellGroupHelper.cpp for detailed rationale.
template<typename Scalar, typename IndexTraits>
template <class AverageRegionalPressureType>
void WellGroupHelper<Scalar, IndexTraits>::
setRegionAveragePressureCalculator(
    const Group& group,
    const FieldPropsManager& fp,
    std::map<std::string, std::unique_ptr<AverageRegionalPressureType>>& regional_average_pressure_calculator)
{
    for (const std::string& groupName : group.groups()) {
        this->setRegionAveragePressureCalculator(this->schedule_.getGroup(groupName, this->report_step_), fp, regional_average_pressure_calculator);
    }
    const auto& gpm = group.gpmaint();
    if (!gpm)
        return;

    const auto& reg = gpm->region();
    if (!reg)
        return;

    if (regional_average_pressure_calculator.count(reg->first) == 0) {
        const std::string name = (reg->first.rfind("FIP", 0) == 0) ? reg->first : "FIP" + reg->first;
        const auto& fipnum = fp.get_int(name);
        regional_average_pressure_calculator[reg->first] =
            std::make_unique<AverageRegionalPressureType>(fipnum);
    }
}

// NOTE: This template member function is defined in the header because the
//       RegionalValues type parameter depends on derived class types that are
//       not available when WellGroupHelper.cpp is compiled. Template functions
//       must be visible at their instantiation point.
//       See WellGroupHelper.cpp for detailed rationale.
template<typename Scalar, typename IndexTraits>
template <class RegionalValues>
void WellGroupHelper<Scalar, IndexTraits>::
updateGpMaintTargetForGroups(const Group& group,
                             const RegionalValues& regional_values,
                             const double dt)
{
    OPM_TIMEFUNCTION();
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_);
        this->updateGpMaintTargetForGroups(group_tmp, regional_values, dt);
    }
    const auto& gpm = group.gpmaint();
    if (!gpm)
        return;

    const auto& region = gpm->region();
    if (!region)
        return;

    const auto [name, number] = *region;
    const Scalar error = gpm->pressure_target() - regional_values.at(name)->pressure(number);
    Scalar current_rate = 0.0;
    const auto& pu = this->phase_usage_info_;
    bool injection = true;
    Scalar sign = 1.0;
    switch (gpm->flow_target()) {
        case GPMaint::FlowTarget::RESV_PROD:
        {
            current_rate = -this->groupState().injection_vrep_rate(group.name());
            injection = false;
            sign = -1.0;
            break;
        }
        case GPMaint::FlowTarget::RESV_OINJ:
        {
            if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
                const auto io = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
                current_rate = this->groupState().injection_reservoir_rates(group.name())[io];
            }
            break;
        }
        case GPMaint::FlowTarget::RESV_WINJ:
        {
            if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
                const auto iw = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
                current_rate = this->groupState().injection_reservoir_rates(group.name())[iw];
            }
            break;
        }
        case GPMaint::FlowTarget::RESV_GINJ:
        {
            if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
                const auto ig = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
                current_rate = this->groupState().injection_reservoir_rates(group.name())[ig];
            }
            break;
        }
        case GPMaint::FlowTarget::SURF_OINJ:
        {
            if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
                const auto io = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
                current_rate = this->groupState().injection_surface_rates(group.name())[io];
            }
            break;
        }
        case GPMaint::FlowTarget::SURF_WINJ:
        {
            if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
                const auto iw = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
                current_rate = this->groupState().injection_surface_rates(group.name())[iw];
            }
            break;
        }
        case GPMaint::FlowTarget::SURF_GINJ:
        {
            if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
                const auto ig = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
                current_rate = this->groupState().injection_surface_rates(group.name())[ig];
            }
            break;
        }
        default:
            throw std::invalid_argument("Invalid Flow target type in GPMAINT");
    }
    auto& gpmaint_state = this->groupState().gpmaint(group.name());
    // we only activate gpmaint if pressure is lower than the target regional pressure for injectors
    // (i.e. error > 0) and higher for producers.
    bool activate = (injection && error > 0) || (!injection && error < 0);
    Scalar rate = 0.0;
    if (activate) {
        rate = gpm->rate(gpmaint_state, current_rate, error, dt);
    } else {
        gpm->resetState(gpmaint_state);
    }
    this->groupState().update_gpmaint_target(group.name(), std::max(Scalar{0.0}, sign * rate));
}

} // namespace Opm

#endif
