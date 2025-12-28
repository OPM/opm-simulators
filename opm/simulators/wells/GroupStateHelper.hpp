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
#ifndef OPM_GROUPSTATEHELPER_HEADER_INCLUDED
#define OPM_GROUPSTATEHELPER_HEADER_INCLUDED

#include <opm/simulators/wells/rescoup/RescoupProxy.hpp>

#include <opm/common/TimingMacros.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/Schedule/Group/GPMaint.hpp>
#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <algorithm>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm
{

template <typename Scalar, typename IndexTraits>
class GroupStateHelper
{
public:
    // RAII guard for temporarily setting wellstate pointer
    class WellStateGuard
    {
    public:
        WellStateGuard(GroupStateHelper& groupStateHelper, WellState<Scalar, IndexTraits>& well_state)
            : groupStateHelper_ {groupStateHelper}
            , previous_state_ptr_ {groupStateHelper_.well_state_}
        {
            // Set the new state directly
            groupStateHelper_.well_state_ = &well_state;
        }

        ~WellStateGuard()
        {
            // Restore the previous state
            groupStateHelper_.well_state_ = previous_state_ptr_;
        }

        // Delete copy and move operations
        WellStateGuard(const WellStateGuard&) = delete;
        WellStateGuard& operator=(const WellStateGuard&) = delete;
        WellStateGuard(WellStateGuard&&) = delete;
        WellStateGuard& operator=(WellStateGuard&&) = delete;

    private:
        GroupStateHelper& groupStateHelper_;
        const WellState<Scalar, IndexTraits>* previous_state_ptr_;
    };

    // RAII guard for temporarily setting groupstate pointer
    class GroupStateGuard
    {
    public:
        GroupStateGuard(GroupStateHelper& group_state_helper, GroupState<Scalar>& group_state)
            : group_state_helper_ {group_state_helper}
            , previous_state_ptr_ {group_state_helper.group_state_}
        {
            // Set the new state directly
            group_state_helper_.group_state_ = &group_state;
        }

        ~GroupStateGuard()
        {
            // Restore the previous state
            group_state_helper_.group_state_ = previous_state_ptr_;
        }

        // Delete copy and move operations
        GroupStateGuard(const GroupStateGuard&) = delete;
        GroupStateGuard& operator=(const GroupStateGuard&) = delete;
        GroupStateGuard(GroupStateGuard&&) = delete;
        GroupStateGuard& operator=(GroupStateGuard&&) = delete;

    private:
        GroupStateHelper& group_state_helper_;
        GroupState<Scalar>* previous_state_ptr_ {nullptr};
    };

    /// @brief RAII guard for scoped DeferredLogger binding
    ///
    /// @details This class ensures that a DeferredLogger pointer is properly set and cleared
    /// in the GroupStateHelper. It follows the RAII pattern to prevent dangling pointer issues
    /// when local DeferredLogger objects go out of scope.
    ///
    /// Supports move semantics to allow returning from factory functions and use with std::optional.
    ///
    /// Lifecycle guarantees:
    /// - Constructor: Sets the DeferredLogger pointer
    /// - Destructor: Clears the pointer (only if not moved-from)
    /// - Move: Transfers ownership; moved-from object becomes inactive
    class ScopedDeferredLoggerGuard
    {
    public:
        ScopedDeferredLoggerGuard(const GroupStateHelper& group_state_helper, DeferredLogger& deferred_logger)
            : group_state_helper_(&group_state_helper)
        {
            group_state_helper_->deferred_logger_ = &deferred_logger;
        }

        ~ScopedDeferredLoggerGuard()
        {
            if (group_state_helper_) {
                group_state_helper_->deferred_logger_ = nullptr;
            }
        }

        // Delete copy operations
        ScopedDeferredLoggerGuard(const ScopedDeferredLoggerGuard&) = delete;
        ScopedDeferredLoggerGuard& operator=(const ScopedDeferredLoggerGuard&) = delete;

        // Enable move operations for std::optional and factory function returns
        ScopedDeferredLoggerGuard(ScopedDeferredLoggerGuard&& other) noexcept
            : group_state_helper_(other.group_state_helper_)
        {
            other.group_state_helper_ = nullptr;
        }

        ScopedDeferredLoggerGuard& operator=(ScopedDeferredLoggerGuard&& other) noexcept
        {
            if (this != &other) {
                if (group_state_helper_) {
                    group_state_helper_->deferred_logger_ = nullptr;
                }
                group_state_helper_ = other.group_state_helper_;
                other.group_state_helper_ = nullptr;
            }
            return *this;
        }

    private:
        const GroupStateHelper* group_state_helper_{nullptr};
    };

    GroupStateHelper(WellState<Scalar, IndexTraits>& well_state,
                    GroupState<Scalar>& group_state,
                    const Schedule& schedule,
                    const SummaryState& summary_state,
                    const GuideRate& guide_rate,
                    const PhaseUsageInfo<IndexTraits>& phase_usage_info,
                    const Parallel::Communication& comm);

    void accumulateGroupEfficiencyFactor(const Group& group, Scalar& factor) const;

    std::pair<bool, Scalar> checkGroupConstraintsInj(const std::string& name,
                                                     const std::string& parent,
                                                     const Group& group,
                                                     const Scalar* rates,
                                                     const Phase injection_phase,
                                                     const Scalar efficiency_factor,
                                                     const std::vector<Scalar>& resv_coeff,
                                                     const bool check_guide_rate) const;

    std::pair<bool, Scalar> checkGroupConstraintsProd(const std::string& name,
                                                      const std::string& parent,
                                                      const Group& group,
                                                      const Scalar* rates,
                                                      const Scalar efficiency_factor,
                                                      const std::vector<Scalar>& resv_coeff,
                                                      const bool check_guide_rate) const;

    const Parallel::Communication& comm() const { return this->comm_; }

    /// @brief Get the deferred logger
    /// @throws assertion failure if no logger has been set via setupScopedDeferredLogger()
    DeferredLogger& deferredLogger() const
    {
        assert(this->deferred_logger_ != nullptr);
        return *this->deferred_logger_;
    }

    std::vector<Scalar> getGroupRatesAvailableForHigherLevelControl(const Group& group, const bool is_injector) const;

    Scalar getGuideRate(const std::string& name, const GuideRateModel::Target target) const;

    Scalar getInjectionGroupTarget(const Group& group,
                                   const Phase& injection_phase,
                                   const std::vector<Scalar>& resv_coeff) const;

    Scalar getProductionGroupTarget(const Group& group, DeferredLogger& deferred_logger) const;

    GuideRate::RateVector getProductionGroupRateVector(const std::string& group_name) const;

    std::optional<Scalar> getWellGroupTargetInjector(const std::string& name,
                                                     const std::string& parent,
                                                     const Group& group,
                                                     const Scalar* rates,
                                                     const Phase injection_phase,
                                                     const Scalar efficiency_factor,
                                                     const std::vector<Scalar>& resv_coeff,
                                                     DeferredLogger& deferred_logger) const;

    std::optional<Scalar> getWellGroupTargetProducer(const std::string& name,
                                                     const std::string& parent,
                                                     const Group& group,
                                                     const Scalar* rates,
                                                     const Scalar efficiency_factor,
                                                     const std::vector<Scalar>& resv_coeff,
                                                     DeferredLogger& deferred_logger) const;

    GuideRate::RateVector getWellRateVector(const std::string& name) const;

    std::vector<std::string> groupChainTopBot(const std::string& bottom, const std::string& top) const;

    /// returns the number of wells that are actively under group control for a given group with name given
    /// by group_name
    int groupControlledWells(const std::string& group_name,
                             const std::string& always_included_child,
                             const bool is_production_group,
                             const Phase injection_phase) const;

    GroupState<Scalar>& groupState() const
    {
        return *this->group_state_;
    }

    const GuideRate& guideRate() const
    {
        return this->guide_rate_;
    }

    bool isRank0() const
    {
        return this->well_state_->isRank0();
    }

    bool isReservoirCouplingMaster() const { return rescoup_.isMaster(); }

    bool isReservoirCouplingMasterGroup(const Group& group) const { return rescoup_.isMasterGroup(group.name()); }

    bool isReservoirCouplingSlave() const { return rescoup_.isSlave(); }

    constexpr int numPhases() const {
        return this->wellState().numPhases();
    }

    int phaseToActivePhaseIdx(const Phase phase) const;

    const PhaseUsageInfo<IndexTraits>& phaseUsage() const {
        return this->phase_usage_info_;
    }

    GroupStateGuard pushGroupState(GroupState<Scalar>& group_state)
    {
        return GroupStateGuard(*this, group_state);
    }

    WellStateGuard pushWellState(WellState<Scalar, IndexTraits>& well_state)
    {
        return WellStateGuard(*this, well_state);
    }

    int reportStepIdx() const
    {
        return report_step_;
    }

    const Schedule& schedule() const
    {
        return this->schedule_;
    }

    ReservoirCoupling::Proxy<Scalar>& rescoup() {
        return rescoup_;
    }

    const ReservoirCoupling::Proxy<Scalar>& rescoup() const {
        return rescoup_;
    }

    ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster() {
        return rescoup_.master();
    }

    const ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster() const {
        return rescoup_.master();
    }

    ReservoirCouplingSlave<Scalar>& reservoirCouplingSlave() {
        return rescoup_.slave();
    }
    const ReservoirCouplingSlave<Scalar>& reservoirCouplingSlave() const {
        return rescoup_.slave();
    }

#ifdef RESERVOIR_COUPLING_ENABLED
    void setReservoirCouplingMaster(ReservoirCouplingMaster<Scalar>* master) {
        rescoup_.setMaster(master);
    }
    void setReservoirCouplingSlave(ReservoirCouplingSlave<Scalar>* slave) {
        rescoup_.setSlave(slave);
    }
#endif

    void setCmodeGroup(const Group& group);

    template <class AverageRegionalPressureType>
    void
    setRegionAveragePressureCalculator(const Group& group,
                                       const FieldPropsManager& fp,
                                       std::map<std::string, std::unique_ptr<AverageRegionalPressureType>>&
                                           regional_average_pressure_calculator) const;

    void setReportStep(int report_step)
    {
        report_step_ = report_step;
    }

    /// @brief Create a scoped guard that binds a DeferredLogger to this helper
    /// @param deferred_logger The logger to bind (must outlive the returned guard)
    /// @return RAII guard that clears the logger on destruction
    ScopedDeferredLoggerGuard setupScopedDeferredLogger(DeferredLogger& deferred_logger) const
    {
        return ScopedDeferredLoggerGuard(*this, deferred_logger);
    }

    const SummaryState& summaryState() const
    {
        return this->summary_state_;
    }

    Scalar sumSolventRates(const Group& group, const bool is_injector) const;

    Scalar sumWellResRates(const Group& group, const int phase_pos, const bool injector) const;

    Scalar sumWellSurfaceRates(const Group& group, const int phase_pos, const bool injector) const;

    Scalar sumWellPhaseRates(bool res_rates,
                             const Group& group,
                             const int phase_pos,
                             const bool injector,
                             const bool network = false) const;

    template <class RegionalValues>
    void updateGpMaintTargetForGroups(const Group& group,
                                      const RegionalValues& regional_values,
                                      const double dt);

    /// update the number of wells that are actively under group control for a given group with name given by
    /// group_name its main usage is to detect cases where there is no wells under group control
    int updateGroupControlledWells(const bool is_production_group,
                                   const Phase injection_phase,
                                   DeferredLogger& deferred_logger);

    void updateGroupProductionRates(const Group& group);

    void updateGroupTargetReduction(const Group& group,
                                    const bool is_injector);

    void updateNetworkLeafNodeProductionRates();

    void updateREINForGroups(const Group& group, bool sum_rank);

    void updateReservoirRatesInjectionGroups(const Group& group);

    void updateState(WellState<Scalar, IndexTraits>& well_state, GroupState<Scalar>& group_state);

    void updateSurfaceRatesInjectionGroups(const Group& group);

    void updateVREPForGroups(const Group& group);

    void updateWellRates(const Group& group,
                         const WellState<Scalar, IndexTraits>& well_state_nupcol,
                         WellState<Scalar, IndexTraits>& well_state) const;

    const WellState<Scalar, IndexTraits>& wellState() const
    {
        return *this->well_state_;
    }

    void updateWellRatesFromGroupTargetScale(const Scalar scale,
                                             const Group& group,
                                             bool is_injector,
                                             WellState<Scalar, IndexTraits>& well_state) const;

    /// Returns the name of the worst offending well and its fraction (i.e. violated_phase / preferred_phase)
    std::pair<std::optional<std::string>, Scalar>
    worstOffendingWell(const Group& group,
                       const Group::ProductionCMode& offended_control,
                       DeferredLogger& deferred_logger) const;

private:
#ifdef RESERVOIR_COUPLING_ENABLED
    /// @brief Convert active phase index to ReservoirCoupling::Phase enum
    /// @param phase_pos Active phase index (0, 1, or 2 in a 3-phase model)
    /// @return The corresponding ReservoirCoupling::Phase enum value
    /// @note This uses the canonical phase ordering (Oil=0, Gas=1, Water=2)
    ReservoirCoupling::Phase activePhaseIdxToRescoupPhase_(int phase_pos) const;
#endif

    //! \brief Compute partial efficiency factor for addback calculation.
    //!
    //! The addback in constraint checking must use the partial efficiency factor
    //! (from local_reduction_level down to the entity), not the accumulated efficiency
    //! (from the entity to the control group). This helper multiplies efficiency
    //! factors from local_reduction_level+1 to the entity.
    //!
    //! \param chain The group chain from control group (top) to entity (bottom)
    //! \param local_reduction_level The level at which addback is applied
    //! \return The partial efficiency factor for addback
    Scalar computeAddbackEfficiency_(const std::vector<std::string>& chain,
                                     std::size_t local_reduction_level) const;

    std::string controlGroup_(const Group& group) const;

    GuideRate::RateVector getGuideRateVector_(const std::vector<Scalar>& rates) const;

    Scalar getReservoirCouplingMasterGroupRate_(const Group& group,
                                                const int phase_pos,
                                                const bool res_rates,
                                                const bool is_injector) const;

    Scalar getSatelliteRate_(const Group& group,
        const int phase_pos,
        const bool res_rates,
        const bool is_injector) const;

    /// check if well/group bottom is a sub well/group of the group top
    bool isInGroupChainTopBot_(const std::string& bottom, const std::string& top) const;

    bool isSatelliteGroup_(const Group& group) const;

    Scalar satelliteInjectionRate_(const ScheduleState& sched,
                                   const Group& group,
                                   const int phase_pos,
                                   bool res_rates) const;

    Scalar satelliteProductionRate_(const ScheduleState& sched,
                                    const Group& group,
                                    const GSatProd::GSatProdGroupProp::Rate rate_comp,
                                    bool res_rates) const;

    std::optional<GSatProd::GSatProdGroupProp::Rate> selectRateComponent_(const int phase_pos) const;

    int updateGroupControlledWellsRecursive_(const std::string& group_name,
                                             const bool is_production_group,
                                             const Phase injection_phase,
                                             DeferredLogger& deferred_logger);

    void updateGroupTargetReductionRecursive_(const Group& group,
                                              const bool is_injector,
                                              std::vector<Scalar>& group_target_reduction);

    const WellState<Scalar, IndexTraits>* well_state_ {nullptr};
    GroupState<Scalar>* group_state_ {nullptr};
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    const GuideRate& guide_rate_;
    // NOTE: The deferred logger does not change the object "meaningful" state, so it should be ok to
    //   make it mutable and store a pointer to it here.
    mutable DeferredLogger* deferred_logger_ {nullptr};
    // NOTE: The phase usage info seems to be read-only throughout the simulation, so it should be safe
    // to store a reference to it here.
    const PhaseUsageInfo<IndexTraits>& phase_usage_info_;
    const Parallel::Communication& comm_;
    int report_step_ {0};
    ReservoirCoupling::Proxy<Scalar> rescoup_{};
};

// -----------------------------------------------------------------------------
// Template member function implementations
// -----------------------------------------------------------------------------

// NOTE: This template member function is defined here in the header because the
//       AverageRegionalPressureType type parameter depends on derived class types that are
//       not available when GroupStateHelper.cpp is compiled. Template functions
//       must be visible at their instantiation point.
//       See GroupStateHelper.cpp for detailed rationale.
template <typename Scalar, typename IndexTraits>
template <class AverageRegionalPressureType>
void
GroupStateHelper<Scalar, IndexTraits>::setRegionAveragePressureCalculator(
    const Group& group,
    const FieldPropsManager& fp,
    std::map<std::string, std::unique_ptr<AverageRegionalPressureType>>& regional_average_pressure_calculator)
    const
{
    for (const std::string& groupName : group.groups()) {
        this->setRegionAveragePressureCalculator(this->schedule_.getGroup(groupName, this->report_step_),
                                                 fp,
                                                 regional_average_pressure_calculator);
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
        regional_average_pressure_calculator[reg->first]
            = std::make_unique<AverageRegionalPressureType>(fipnum);
    }
}

// NOTE: This template member function is defined in the header because the
//       RegionalValues type parameter depends on derived class types that are
//       not available when GroupStateHelper.cpp is compiled. Template functions
//       must be visible at their instantiation point.
//       See GroupStateHelper.cpp for detailed rationale.
template <typename Scalar, typename IndexTraits>
template <class RegionalValues>
void
GroupStateHelper<Scalar, IndexTraits>::updateGpMaintTargetForGroups(const Group& group,
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
    case GPMaint::FlowTarget::RESV_PROD: {
        current_rate = -this->groupState().injection_vrep_rate(group.name());
        injection = false;
        sign = -1.0;
        break;
    }
    case GPMaint::FlowTarget::RESV_OINJ: {
        if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
            const auto io = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            current_rate = this->groupState().injection_reservoir_rates(group.name())[io];
        }
        break;
    }
    case GPMaint::FlowTarget::RESV_WINJ: {
        if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            const auto iw = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            current_rate = this->groupState().injection_reservoir_rates(group.name())[iw];
        }
        break;
    }
    case GPMaint::FlowTarget::RESV_GINJ: {
        if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
            const auto ig = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
            current_rate = this->groupState().injection_reservoir_rates(group.name())[ig];
        }
        break;
    }
    case GPMaint::FlowTarget::SURF_OINJ: {
        if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
            const auto io = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            current_rate = this->groupState().injection_surface_rates(group.name())[io];
        }
        break;
    }
    case GPMaint::FlowTarget::SURF_WINJ: {
        if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            const auto iw = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            current_rate = this->groupState().injection_surface_rates(group.name())[iw];
        }
        break;
    }
    case GPMaint::FlowTarget::SURF_GINJ: {
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
    this->groupState().update_gpmaint_target(group.name(), std::max(Scalar {0.0}, sign * rate));
}

} // namespace Opm

#endif  // OPM_GROUPSTATE_HELPER_HEADER_INCLUDED
