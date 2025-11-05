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

#ifndef OPM_GROUP_TARGET_CALCULATOR_HPP
#define OPM_GROUP_TARGET_CALCULATOR_HPP
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/FractionCalculator.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/WellGroupHelper.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <variant>
#include <optional>
#include <vector>
#include <string>

namespace Opm {

/**
 * Calculate group-level targets for production and injection.
 *
 * This class traverses the group hierarchy to determine effective control
 * modes and targets, applying guide-rate based distribution, sales limits
 * (e.g., GCONSALE), efficiency factors, and RESV coefficients where
 * applicable. It provides a uniform interface for both producers and injectors
 * and consolidates common logic through nested helper classes.
 */
template<class Scalar, class IndexTraits>
class GroupTargetCalculator {
public:
    /**
     * Union of control-mode types used by group target calculations.
     * Holds Group::InjectionCMode for injection, Group::ProductionCMode for
     * production; std::monostate denotes that no specific control applies.
     */
    using ControlMode = std::variant<
        std::monostate,
        Group::InjectionCMode,
        Group::ProductionCMode
    >;
    using FractionCalculator = WGHelpers::FractionCalculator<Scalar, IndexTraits>;
    using InjectionTargetCalculator = WGHelpers::InjectionTargetCalculator<Scalar, IndexTraits>;
    using TargetCalculator = WGHelpers::TargetCalculator<Scalar, IndexTraits>;
    using WellGroupHelperType = WellGroupHelper<Scalar, IndexTraits>;

    /** Generic result for a computed target and its control mode. */
    struct TargetInfo {
        Scalar target;
        ControlMode cmode;
    };
    /** Result for a production target with its production control mode. */
    struct ProductionTargetInfo {
        Scalar target;
        Group::ProductionCMode cmode;
    };
    /** Result for an injection target with its injection control mode. */
    struct InjectionTargetInfo {
        Scalar target;
        Group::InjectionCMode cmode;
    };

    /**
     * Shared logic for injector and producer paths.
     *
     * Provides the common context (schedule, states, guide rates, phase usage,
     * PVT/FIP regions, logger) and computes a target for the requested bottom
     * group. The optional injection_phase indicates whether the calculation is
     * for injection (has value) or production (no value). Internally it
     * selects phase-dependent RESV coefficients and performs a recursive
     * traversal to accumulate limits and select the effective control mode.
     */
    class GeneralCalculator {
    public:
        GeneralCalculator(
            GroupTargetCalculator& calculator,
            const Group& original_group,  // the bottom group we want to calculate the target for
            std::optional<ReservoirCoupling::Phase> injection_phase = std::nullopt
        );
        const Group& originalGroup() const { return this->original_group_; }
        std::optional<TargetInfo> calculateGroupTarget();
        DeferredLogger& deferredLogger() { return this->parent_calculator_.deferredLogger(); }
        int fipnum() const { return this->parent_calculator_.fipnum(); }
        const GroupState<Scalar>& groupState() const { return this->parent_calculator_.groupState(); }
        const GuideRate& guideRate() const { return this->parent_calculator_.guideRate(); }
        Phase injectionPhase_() const;
        bool isInjector() const { return this->injection_phase_.has_value(); }
        bool isProducer() const { return !this->injection_phase_.has_value(); }
        const PhaseUsageInfo<IndexTraits>& phaseUsage() const { return this->parent_calculator_.phaseUsage(); }
        int pvtreg() const { return this->parent_calculator_.pvtreg(); }
        int reportStepIdx() const { return this->parent_calculator_.reportStepIdx(); }
        const std::vector<Scalar>& resvCoeffsInj() const { return this->parent_calculator_.resvCoeffsInj(); }
        const std::vector<Scalar>& resvCoeffsProd() const { return this->resv_coeffs_prod_; }
        const Schedule& schedule() const { return this->parent_calculator_.schedule(); }
        const SummaryState& summaryState() const { return this->parent_calculator_.summaryState(); }
        const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel() const {
            return this->parent_calculator_.wellModel();
        }
        const WellState<Scalar, IndexTraits>& wellState() const { return this->parent_calculator_.wellState(); }
        const WellGroupHelper<Scalar, IndexTraits>& wgHelper() const { return this->parent_calculator_.wgHelper(); }
    private:
        std::optional<TargetInfo> calculateGroupTargetRecursive_(const Group& group, const Scalar efficiency_factor);
        bool hasFldOrNoneControl_(const Group& group) const;
        const Group& parentGroup(const Group& group) const {
            return this->schedule().getGroup(group.parent(), this->reportStepIdx());
        }
        bool parentGroupControlAvailable_(const Group& group) const;
        Phase reservoirCouplingToOpmPhase_(ReservoirCoupling::Phase reservoir_coupling_phase) const;

        GroupTargetCalculator& parent_calculator_;
        const Group& original_group_; // The bottom group we want to calculate the target for
        std::optional<ReservoirCoupling::Phase> injection_phase_;
        std::vector<Scalar> resv_coeffs_prod_;
    };

    /**
     * Distribute a top-level target down to the requested group.
     *
     * After the effective top target and control mode are known, this helper
     * walks the chain of groups from the top to the bottom group and applies
     * either a production TargetCalculator or an InjectionTargetCalculator.
     * A FractionCalculator is used where guide-rate fractions/reductions are
     * required.
     */
    class TopToBottomCalculator {
    public:
        TopToBottomCalculator(
            GeneralCalculator& parent_calculator,
            const Group& group,
            Scalar efficiency_factor
        );

        const Group& bottomGroup() const { return this->parent_calculator_.originalGroup(); }
        std::optional<TargetInfo> calculateGroupTarget();
        DeferredLogger& deferredLogger() { return this->parent_calculator_.deferredLogger(); }
        const GConSale& gconsale() const {
            return this->schedule()[this->reportStepIdx()].gconsale();
        }
        const GroupState<Scalar>& groupState() const { return this->parent_calculator_.groupState(); }
        const GuideRate& guideRate() const { return this->parent_calculator_.guideRate(); }
        Phase injectionPhase_() const { return this->parent_calculator_.injectionPhase_(); }
        bool isInjector() const { return this->parent_calculator_.isInjector(); }
        bool isProducer() const { return this->parent_calculator_.isProducer(); }
        const PhaseUsageInfo<IndexTraits>& phaseUsage() const { return this->parent_calculator_.phaseUsage(); }
        int reportStepIdx() const { return this->parent_calculator_.reportStepIdx(); }
        const std::vector<Scalar>& resvCoeffsInj() const { return this->parent_calculator_.resvCoeffsInj(); }
        const std::vector<Scalar>& resvCoeffsProd() const { return this->parent_calculator_.resvCoeffsProd(); }
        const Schedule& schedule() const { return this->parent_calculator_.schedule(); }
        const SummaryState& summaryState() const { return this->parent_calculator_.summaryState(); }
        const WellState<Scalar, IndexTraits>& wellState() const { return this->parent_calculator_.wellState(); }
        const WellGroupHelper<Scalar, IndexTraits>& wgHelper() const { return this->parent_calculator_.wgHelper(); }

    private:
        Scalar getTopLevelTarget_();
        Scalar getGratSalesInjectionTarget_() const;
        Scalar getGratSalesProductionTarget_() const;
        std::vector<std::string> getGroupChainTopBot_() const;
        void initForInjector_();
        void initForProducer_();
        Scalar localFraction_(const std::string& group_name);
        Scalar localReduction_(const std::string& group_name);

        GeneralCalculator& parent_calculator_;
        const Group& group_;
        Scalar efficiency_factor_;
        // Active calculator used for distributing the target along the chain:
        // either production TargetCalculator or InjectionTargetCalculator.
        std::variant<std::monostate, TargetCalculator, InjectionTargetCalculator> target_calculator_;
        // Since FractionCalculator does not have a default constructor, we use std::optional
        // to conditionally initialize it based on whether we are dealing with an injector or producer.
        std::optional<FractionCalculator> fraction_calculator_;
        ControlMode toplevel_control_mode_;
    };

    /**
     * Construct a calculator bound to one report step and simulator state.
     */
    GroupTargetCalculator(
        const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model,
        const WellGroupHelper<Scalar, IndexTraits>& wg_helper,
        DeferredLogger& deferred_logger
    );
    DeferredLogger& deferredLogger() { return this->deferred_logger_; }
    int fipnum() const { return this->fipnum_; }
    /** Compute injection target for group in the given injection phase. */
    std::optional<InjectionTargetInfo> groupInjectionTarget(
        const Group& group, ReservoirCoupling::Phase injection_phase
    );
    /** Compute production target for group. */
    std::optional<ProductionTargetInfo> groupProductionTarget(const Group& group);
    const GroupState<Scalar>& groupState() const { return this->group_state_; }
    const GuideRate& guideRate() const { return this->guide_rate_; }
    const PhaseUsageInfo<IndexTraits>& phaseUsage() const { return this->phase_usage_; }
    int pvtreg() const { return this->pvtreg_; }
    int reportStepIdx() const { return this->report_step_idx_; }
    const std::vector<Scalar>& resvCoeffsInj() const { return this->resv_coeffs_inj_; }
    const Schedule& schedule() const { return this->schedule_; }
    const SummaryState& summaryState() const { return this->summary_state_; }
    const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel() const { return this->well_model_; }
    const WellState<Scalar, IndexTraits>& wellState() const { return this->well_state_; }
    const WellGroupHelper<Scalar, IndexTraits>& wgHelper() const { return this->wg_helper_; }
private:
    const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model_;
    const WellGroupHelper<Scalar, IndexTraits>& wg_helper_;
    const WellState<Scalar, IndexTraits >& well_state_;
    const GroupState<Scalar>& group_state_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    const PhaseUsageInfo<IndexTraits>& phase_usage_;
    const GuideRate& guide_rate_;
    int report_step_idx_;
    std::vector<Scalar> resv_coeff_;
    DeferredLogger& deferred_logger_;
    int fipnum_; // FIP region for the groups
    int pvtreg_; // PVT region for the groups
    std::vector<Scalar> resv_coeffs_inj_;
};

} // namespace Opm

#endif // OPM_GROUP_TARGET_CALCULATOR_HPP
