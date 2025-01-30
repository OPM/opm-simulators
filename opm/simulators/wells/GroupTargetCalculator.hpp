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
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <variant>
#include <optional>
#include <vector>
#include <string>

namespace Opm {

template<class Scalar>
class GroupTargetCalculator {
public:
    using ControlMode = std::variant<
        std::monostate,
        Group::InjectionCMode,
        Group::ProductionCMode
    >;
    using FractionCalculator = WGHelpers::FractionCalculator<Scalar>;
    using InjectionTargetCalculator = WGHelpers::InjectionTargetCalculator<Scalar>;
    using TargetCalculator = WGHelpers::TargetCalculator<Scalar>;
    struct TargetInfo {
        Scalar target;
        ControlMode cmode;
    };
    struct ProductionTargetInfo {
        Scalar target;
        Group::ProductionCMode cmode;
    };
    struct InjectionTargetInfo {
        Scalar target;
        Group::InjectionCMode cmode;
    };

    // Generalizes calculation for both injectors and producers to avoid code duplication.
    class GeneralCalculator {
    public:
        GeneralCalculator(
            GroupTargetCalculator& calculator,
            const Group& original_group,  // the bottom group we want to calculate the target for
            std::optional<Phase> injection_phase = std::nullopt
        );
        const Group& originalGroup() const { return this->original_group_; }
        std::optional<TargetInfo> calculateGroupTarget();
        DeferredLogger& deferredLogger() { return this->parent_calculator_.deferredLogger(); }
        int fipnum() const { return this->parent_calculator_.fipnum(); }
        const GroupState<Scalar>& groupState() const { return this->parent_calculator_.groupState(); }
        const GuideRate& guideRate() const { return this->parent_calculator_.guideRate(); }
        Phase injectionPhase() const { return this->injection_phase_.value(); }
        bool isInjector() const { return this->injection_phase_.has_value(); }
        bool isProducer() const { return !this->injection_phase_.has_value(); }
        const PhaseUsage& phaseUsage() const { return this->parent_calculator_.phaseUsage(); }
        int pvtreg() const { return this->parent_calculator_.pvtreg(); }
        int reportStepIdx() const { return this->parent_calculator_.reportStepIdx(); }
        const std::vector<Scalar>& resvCoeffsInj() const { return this->parent_calculator_.resvCoeffsInj(); }
        const std::vector<Scalar>& resvCoeffsProd() const { return this->resv_coeffs_prod_; }
        const Schedule& schedule() const { return this->parent_calculator_.schedule(); }
        const SummaryState& summaryState() const { return this->parent_calculator_.summaryState(); }
        const BlackoilWellModelGeneric<Scalar>& wellModel() const {
            return this->parent_calculator_.wellModel();
        }
        const WellState<Scalar>& wellState() const { return this->parent_calculator_.wellState(); }
    private:
        std::optional<TargetInfo> calculateGroupTargetRecursive_(const Group& group, const Scalar efficiency_factor);
        bool hasFldOrNoneControl_(const Group& group) const;
        const Group& parentGroup(const Group& group) const {
            return this->schedule().getGroup(group.parent(), this->reportStepIdx());
        }
        bool parentGroupControlAvailable_(const Group& group) const;

        GroupTargetCalculator& parent_calculator_;
        const Group& original_group_; // The bottom group we want to calculate the target for
        std::optional<Phase> injection_phase_;
        std::vector<Scalar> resv_coeffs_prod_;
    };

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
        Phase injectionPhase() const { return this->parent_calculator_.injectionPhase(); }
        bool isInjector() const { return this->parent_calculator_.isInjector(); }
        bool isProducer() const { return this->parent_calculator_.isProducer(); }
        const PhaseUsage& phaseUsage() const { return this->parent_calculator_.phaseUsage(); }
        int reportStepIdx() const { return this->parent_calculator_.reportStepIdx(); }
        const std::vector<Scalar>& resvCoeffsInj() const { return this->parent_calculator_.resvCoeffsInj(); }
        const std::vector<Scalar>& resvCoeffsProd() const { return this->parent_calculator_.resvCoeffsProd(); }
        const Schedule& schedule() const { return this->parent_calculator_.schedule(); }
        const SummaryState& summaryState() const { return this->parent_calculator_.summaryState(); }
        const WellState<Scalar>& wellState() const { return this->parent_calculator_.wellState(); }

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
        std::variant<std::monostate, TargetCalculator, InjectionTargetCalculator> target_calculator_;
        // Since FractionCalculator does not have a default constructor, we use std::optional
        // to conditionally initialize it based on whether we are dealing with an injector or producer.
        std::optional<FractionCalculator> fraction_calculator_;
        ControlMode toplevel_control_mode_;
    };

    GroupTargetCalculator(
        const BlackoilWellModelGeneric<Scalar>& well_model,
        const WellState<Scalar>& well_state,
        const GroupState<Scalar>& group_state,
        const Schedule& schedule,
        const SummaryState& summary_state,
        const PhaseUsage& phase_usage,
        const GuideRate& guide_rate,
        const int report_step_idx,
        DeferredLogger& deferred_logger
    );
    DeferredLogger& deferredLogger() { return this->deferred_logger_; }
    int fipnum() const { return this->fipnum_; }
    std::optional<InjectionTargetInfo> groupInjectionTarget(const Group& group, Phase injection_phase);
    std::optional<ProductionTargetInfo> groupProductionTarget(const Group& group);
    const GroupState<Scalar>& groupState() const { return this->group_state_; }
    const GuideRate& guideRate() const { return this->guide_rate_; }
    const PhaseUsage& phaseUsage() const { return this->phase_usage_; }
    int pvtreg() const { return this->pvtreg_; }
    int reportStepIdx() const { return this->report_step_idx_; }
    const std::vector<Scalar>& resvCoeffsInj() const { return this->resv_coeffs_inj_; }
    const Schedule& schedule() const { return this->schedule_; }
    const SummaryState& summaryState() const { return this->summary_state_; }
    const BlackoilWellModelGeneric<Scalar>& wellModel() const { return this->well_model_; }
    const WellState<Scalar>& wellState() const { return this->well_state_; }

private:
    const BlackoilWellModelGeneric<Scalar>& well_model_;
    const WellState<Scalar>& well_state_;
    const GroupState<Scalar>& group_state_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    const PhaseUsage& phase_usage_;
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
