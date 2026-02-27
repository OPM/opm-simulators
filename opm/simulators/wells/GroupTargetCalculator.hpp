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
#include <opm/common/ErrorMacros.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/FractionCalculator.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
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
    using FractionCalculator = GroupStateHelpers::FractionCalculator<Scalar, IndexTraits>;
    using InjectionTargetCalculator = GroupStateHelpers::InjectionTargetCalculator<Scalar, IndexTraits>;
    using TargetCalculator = GroupStateHelpers::TargetCalculator<Scalar, IndexTraits>;
    using GroupStateHelperType = GroupStateHelper<Scalar, IndexTraits>;

    enum class TargetType {
        Injection,
        Production
    };

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
    /** Per-rate-type production constraints for a group.
     *
     *  A master group defines both an active control-mode target and limits for
     *  other rate types (ORAT, WRAT, GRAT, LRAT, RESV).  The active target
     *  alone is not sufficient for reservoir-coupling slaves: a slave group must
     *  know every effective limit so that it can enforce the most restrictive
     *  constraint for each rate type independently (the active cmode on the
     *  master side may differ from what is binding on the slave side).
     *
     *  Each limit field holds the guide-rate-distributed value that applies
     *  to this group, or -1 when no ancestor in the hierarchy defines a
     *  limit for that rate type.
     */
    struct ProductionConstraintsInfo {
        Scalar active_target;
        Group::ProductionCMode active_cmode;
        Scalar oil_limit = -1;
        Scalar water_limit = -1;
        Scalar gas_limit = -1;
        Scalar liquid_limit = -1;
        Scalar resv_limit = -1;
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
        using TargetCalculatorType = std::variant<std::monostate, TargetCalculator, InjectionTargetCalculator>;
        GeneralCalculator(
            GroupTargetCalculator& calculator,
            const Group& original_group,  // the bottom group we want to calculate the target for
            std::optional<ReservoirCoupling::Phase> injection_phase = std::nullopt
        );
        /// Construct for computing a limit for an explicit production rate type.
        /// The recursion stops at ancestors with has_control(explicit_cmode).
        GeneralCalculator(
            GroupTargetCalculator& calculator,
            const Group& original_group,
            Group::ProductionCMode explicit_cmode
        );
        std::optional<TargetInfo> calculateGroupTarget();
        DeferredLogger& deferredLogger() { return this->parent_calculator_.deferredLogger(); }
        // Const overload: allows logging from const methods.
        DeferredLogger& deferredLogger() const { return this->parent_calculator_.deferredLogger(); }
        int fipnum() const { return this->parent_calculator_.fipnum(); }
        const GConSale& gconsale() const {
            return this->schedule()[this->reportStepIdx()].gconsale();
        }
        const GroupState<Scalar>& groupState() const { return this->parent_calculator_.groupState(); }
        TargetCalculatorType getInjectionTargetCalculator(const Group& group);
        TargetCalculatorType getProductionTargetCalculator(const Group& group) const;
        TargetCalculatorType getTargetCalculator(const Group& group);
        TargetInfo getGroupTargetNoGuideRate(const Group& group);
        const GuideRate& guideRate() const { return this->parent_calculator_.guideRate(); }
        bool hasGuideRate(const Group& group) const { return this->hasGuideRate(group.name()); }
        bool hasGuideRate(const std::string& name) const {
            if (this->targetType() == TargetType::Injection) {
                return this->guideRate().has(name, this->injectionPhase_());
            }
            return this->guideRate().has(name);
        }
        bool hasFldOrNoneControl(const Group& group);
        Phase injectionPhase_() const;
        const Group& originalGroup() const { return this->original_group_; }
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
        TargetType targetType() const {
            return this->injection_phase_.has_value() ? TargetType::Injection : TargetType::Production;
        }
        const WellState<Scalar, IndexTraits>& wellState() const { return this->parent_calculator_.wellState(); }
        const GroupStateHelperType& groupStateHelper() const { return this->parent_calculator_.groupStateHelper(); }
        /// Returns the explicit production cmode if set (for per-rate-type limit computation).
        std::optional<Group::ProductionCMode> explicitProdCmode() const { return this->explicit_prod_cmode_; }
    private:
        std::optional<TargetInfo> calculateGroupTargetRecursive_(const Group& group, const Scalar efficiency_factor);
        const Group& parentGroup(const Group& group) const {
            return this->schedule().getGroup(group.parent(), this->reportStepIdx());
        }
        bool parentGroupControlAvailable_(const Group& group);
        Phase reservoirCouplingToOpmPhase_(ReservoirCoupling::Phase reservoir_coupling_phase) const;

        GroupTargetCalculator& parent_calculator_;
        const Group& original_group_; // The bottom group we want to calculate the target for
        std::optional<ReservoirCoupling::Phase> injection_phase_;
        std::optional<Group::ProductionCMode> explicit_prod_cmode_;
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
        using TargetCalculatorType = std::variant<std::monostate, TargetCalculator, InjectionTargetCalculator>;

        constexpr static Scalar TARGET_RATE_TOLERANCE = 1e-12;

        TopToBottomCalculator(
            GeneralCalculator& parent_calculator,
            const Group& top_group,
            const Group& bottom_group,
            Scalar efficiency_factor
        );

        std::optional<TargetInfo> calculateGroupTarget();
        DeferredLogger& deferredLogger() { return this->parent_calculator_.deferredLogger(); }
        // Const overload: allows logging from const methods (logical constness for external logger).
        DeferredLogger& deferredLogger() const { return this->parent_calculator_.deferredLogger(); }
        const GroupState<Scalar>& groupState() const { return this->parent_calculator_.groupState(); }
        const GuideRate& guideRate() const { return this->parent_calculator_.guideRate(); }
        TargetType targetType() const { return this->parent_calculator_.targetType(); }
        const PhaseUsageInfo<IndexTraits>& phaseUsage() const { return this->parent_calculator_.phaseUsage(); }
        int reportStepIdx() const { return this->parent_calculator_.reportStepIdx(); }
        const std::vector<Scalar>& resvCoeffsInj() const { return this->parent_calculator_.resvCoeffsInj(); }
        const std::vector<Scalar>& resvCoeffsProd() const { return this->parent_calculator_.resvCoeffsProd(); }
        const Schedule& schedule() const { return this->parent_calculator_.schedule(); }
        const SummaryState& summaryState() const { return this->parent_calculator_.summaryState(); }
        const WellState<Scalar, IndexTraits>& wellState() const { return this->parent_calculator_.wellState(); }
        const GroupStateHelperType& groupStateHelper() const { return this->parent_calculator_.groupStateHelper(); }

    private:
        bool bottomGroupHasIndividualControl_();
        Scalar computeAddbackEfficiency_(const std::vector<std::string>& chain,
                                         const std::size_t local_reduction_level) const;
        Scalar getBottomGroupCurrentRateAvailable_() const;
        std::vector<std::string> getGroupChainTopBot_() const;
        std::size_t getLocalReductionLevel_(const std::vector<std::string>& chain);
        TargetCalculatorType getProductionTargetCalculator_(const Group& group) const {
            return this->parent_calculator_.getProductionTargetCalculator(group); }
        TargetCalculatorType getInjectionTargetCalculator_(const Group& group) const {
            return this->parent_calculator_.getInjectionTargetCalculator(group); }
        TargetCalculatorType getInjectionTargetCalculator(const Group& group) const;
        TargetCalculatorType getProductionTargetCalculator(const Group& group) const;
        /// @brief Get the slave group's total reservoir rate for RESV scaling
        /// @details For reservoir coupling master groups, this returns the sum of all phase
        ///          reservoir rates from the slave (in slave's reservoir units). This is used
        ///          when scaling targets in RESV mode.
        /// @param group The master group to get the corresponding slave group reservoir rate for
        /// @return Total reservoir rate, or throws an error if not a RC master group
        Scalar getSlaveGroupReservoirRate_(const Group& master_group);
        TargetInfo getGroupTargetNoGuideRate_(const Group& group) const {
            return this->parent_calculator_.getGroupTargetNoGuideRate(group);
        }
        ControlMode getToplevelControlMode_() const;
        Scalar getTopLevelTarget_();
        bool hasFldOrNoneControl_(const Group& group) {
            return this->parent_calculator_.hasFldOrNoneControl(group);
        }
        bool hasFLDControl_(const Group& group) const;
        bool hasGuideRate_(const std::string& name) const {
            return this->parent_calculator_.hasGuideRate(name);
        }
        void initForInjector_();
        void initForProducer_();
        Phase injectionPhase_() const { return this->parent_calculator_.injectionPhase_(); }
        bool isProducerAndRESVControl_(const Group& group) const;
        Scalar localFraction_(const std::string& group_name);
        Scalar localReduction_(const std::string& group_name);

        GeneralCalculator& parent_calculator_;
        const Group& top_group_;
        const Group& bottom_group_;
        // Accumulated efficiency factor along the chain from top to bottom, excluding the top group.
        Scalar chain_efficiency_factor_;
        // Active calculator used for distributing the target along the chain:
        // either production TargetCalculator or InjectionTargetCalculator.
        std::variant<std::monostate, TargetCalculator, InjectionTargetCalculator> target_calculator_;
        // Since FractionCalculator does not have a default constructor, we use std::optional
        // to conditionally initialize it based on whether we are dealing with an injector or producer.
        std::optional<FractionCalculator> fraction_calculator_;
    };

    /**
     * Construct a calculator bound to one report step and simulator state.
     */
    GroupTargetCalculator(
        const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model,
        const GroupStateHelperType& group_state_helper
    );
    DeferredLogger& deferredLogger() { return this->group_state_helper_.deferredLogger(); }
    // Const overload: allows logging from const methods.
    DeferredLogger& deferredLogger() const { return this->group_state_helper_.deferredLogger(); }
    int fipnum() const { return this->fipnum_; }
    /** Compute injection target for group in the given injection phase. */
    std::optional<InjectionTargetInfo> groupInjectionTarget(
        const Group& group, ReservoirCoupling::Phase injection_phase
    );
    /** Compute production target for group. */
    std::optional<ProductionTargetInfo> groupProductionTarget(const Group& group);
    /** Compute per-rate-type production constraints for a group.
     *  Returns the active target plus effective limits for all rate types
     *  with defined limits in the group hierarchy. */
    std::optional<ProductionConstraintsInfo> groupProductionConstraints(const Group& group);
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
    const GroupStateHelperType& groupStateHelper() const { return this->group_state_helper_; }
private:
    const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model_;
    const GroupStateHelperType& group_state_helper_;
    const WellState<Scalar, IndexTraits >& well_state_;
    const GroupState<Scalar>& group_state_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    const PhaseUsageInfo<IndexTraits>& phase_usage_;
    const GuideRate& guide_rate_;
    int report_step_idx_;
    std::vector<Scalar> resv_coeff_;
    int fipnum_; // FIP region for the groups
    int pvtreg_; // PVT region for the groups
    std::vector<Scalar> resv_coeffs_inj_;
};

} // namespace Opm

#endif // OPM_GROUP_TARGET_CALCULATOR_HPP
