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
#include <opm/input/eclipse/Schedule/ResCoup/GrupSlav.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/Schedule/Group/GPMaint.hpp>
#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/gatherDeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
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

    /// @brief RAII guard that owns a DeferredLogger and auto-gathers on destruction
    ///
    /// @details This class provides a complete lifecycle for deferred logging in parallel
    /// simulations. It owns a DeferredLogger instance and automatically handles:
    /// - Pushing the logger onto the stack (saving the previous logger)
    /// - On destruction: gathering messages across MPI ranks, logging on rank 0
    ///   (if terminal_output enabled), and restoring the previous logger
    ///
    /// Usage:
    ///   auto guard = groupStateHelper.pushLogger();
    ///   // Use groupStateHelper.deferredLogger() to log messages
    ///   // On scope exit: gather, log (if terminal), restore previous logger
    class ScopedLoggerGuard
    {
    public:
        /// @brief Constructor for scoped logger guard
        /// @param helper Reference to the GroupStateHelper
        /// @param do_mpi_gather If true, gather messages across MPI ranks on destruction.
        ///        Set to false when called from contexts where MPI synchronization is not
        ///        possible (e.g., NLDD domain-local operations).
        explicit ScopedLoggerGuard(const GroupStateHelper& helper, bool do_mpi_gather = true)
            : helper_(&helper)
            , previous_(helper.deferred_logger_)
            , do_mpi_gather_(do_mpi_gather)
        {
            helper_->deferred_logger_ = &logger_;
        }

        ~ScopedLoggerGuard()
        {
            if (helper_) {
                if (do_mpi_gather_) {
                    // 1. Gather messages across MPI ranks
                    DeferredLogger global = gatherDeferredLogger(logger_, helper_->comm());

                    // 2. Log on rank 0 (if terminal_output enabled)
                    if (helper_->terminalOutput()) {
                        global.logMessages();
                    }
                } else {
                    // Just log locally without MPI gather
                    if (helper_->terminalOutput()) {
                        logger_.logMessages();
                    }
                }

                // 3. Restore previous logger
                helper_->deferred_logger_ = previous_;
            }
        }

        // Delete copy operations and move assignment
        ScopedLoggerGuard(const ScopedLoggerGuard&) = delete;
        ScopedLoggerGuard& operator=(const ScopedLoggerGuard&) = delete;
        ScopedLoggerGuard& operator=(ScopedLoggerGuard&&) = delete;

        // Move constructor required for pushLogger() return-by-value (must be
        // available even if elided by RVO)
        ScopedLoggerGuard(ScopedLoggerGuard&& other) noexcept
            : helper_(other.helper_)
            , logger_(std::move(other.logger_))
            , previous_(other.previous_)
            , do_mpi_gather_(other.do_mpi_gather_)
        {
            // Update the helper's pointer to our moved logger
            if (helper_) {
                helper_->deferred_logger_ = &logger_;
            }
            other.helper_ = nullptr;
            other.previous_ = nullptr;
        }

    private:
        // Pointer (not reference) to allow nulling in move constructor
        const GroupStateHelper* helper_{nullptr};
        DeferredLogger logger_;            // Owned logger
        DeferredLogger* previous_{nullptr}; // For restore
        bool do_mpi_gather_{true};         // Whether to gather messages across MPI ranks
    };

    using GroupTarget = typename SingleWellState<Scalar, IndexTraits>::GroupTarget;

    GroupStateHelper(WellState<Scalar, IndexTraits>& well_state,
                    GroupState<Scalar>& group_state,
                    const Schedule& schedule,
                    const SummaryState& summary_state,
                    const GuideRate& guide_rate,
                    const PhaseUsageInfo<IndexTraits>& phase_usage_info,
                    const Parallel::Communication& comm,
                    bool terminal_output);

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

    std::pair<Group::ProductionCMode, Scalar>
    checkGroupProductionConstraints(const Group& group) const;

    const Parallel::Communication& comm() const { return this->comm_; }

    /// @brief Get the deferred logger
    /// @throws std::logic_error if no logger has been set via pushLogger()
    DeferredLogger& deferredLogger() const
    {
        if (this->deferred_logger_ == nullptr) {
            throw std::logic_error("DeferredLogger not set. Call pushLogger() first.");
        }
        return *this->deferred_logger_;
    }

    std::vector<Scalar> getGroupRatesAvailableForHigherLevelControl(const Group& group, const bool is_injector) const;

    Scalar getInjectionGroupTarget(const Group& group,
                                   const Phase& injection_phase,
                                   const std::vector<Scalar>& resv_coeff) const;

    /// @brief Get the guide rate target mode for an injection phase
    /// @param injection_phase The injection phase (WATER, OIL, or GAS)
    /// @return The corresponding GuideRateModel::Target for the injection phase
    GuideRateModel::Target getInjectionGuideTargetMode(Phase injection_phase) const;

    Scalar getProductionGroupTarget(const Group& group) const;

    /// Get the production target for a specific control mode (not necessarily the active one).
    Scalar getProductionGroupTargetForMode(const Group& group,
                                           Group::ProductionCMode cmode) const;

    /// @brief Get the guide rate target mode for a production group
    /// @param group The production group
    /// @return The GuideRateModel::Target based on the group's production control mode
    GuideRateModel::Target getProductionGuideTargetMode(const Group& group) const;

    std::pair<Scalar, Group::ProductionCMode>
    getAutoChokeGroupProductionTargetRate(const Group& bottom_group,
                                          const Group& group,
                                          const std::vector<Scalar>& resv_coeff,
                                          Scalar efficiencyFactor) const;

    GuideRate::RateVector getProductionGroupRateVector(const std::string& group_name) const;

    std::optional<GroupTarget> getWellGroupTargetInjector(const std::string& name,
                                                          const std::string& parent,
                                                          const Group& group,
                                                          const Scalar* rates,
                                                          const Phase injection_phase,
                                                          const Scalar efficiency_factor,
                                                          const std::vector<Scalar>& resv_coeff) const;

    std::optional<GroupTarget> getWellGroupTargetProducer(const std::string& name,
                                                          const std::string& parent,
                                                          const Group& group,
                                                          const Scalar* rates,
                                                          const Scalar efficiency_factor,
                                                          const std::vector<Scalar>& resv_coeff) const;

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

    bool isReservoirCouplingSlaveGroup(const Group& group) const { return rescoup_.isSlaveGroup(group.name()); }

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

    /// @brief Push a new logger onto the stack with auto-cleanup on destruction
    ///
    /// @details Creates a new DeferredLogger and pushes it onto the logger stack.
    /// When the returned guard goes out of scope:
    /// 1. If do_mpi_gather is true: Messages are gathered across MPI ranks via gatherDeferredLogger()
    ///    and logged on rank 0 (if terminal_output is enabled)
    /// 2. If do_mpi_gather is false: Messages are logged locally without MPI gather
    /// 3. The previous logger is restored
    ///
    /// @param do_mpi_gather If true (default), gather messages across MPI ranks on destruction.
    ///        Set to false when called from contexts where MPI synchronization is not
    ///        possible (e.g., NLDD domain-local operations called at different times on different ranks).
    /// @return RAII guard that owns the logger and handles cleanup
    ScopedLoggerGuard pushLogger(bool do_mpi_gather = true) const
    {
        return ScopedLoggerGuard(*this, do_mpi_gather);
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

    bool terminalOutput() const
    {
         return this->terminal_output_;
    }

    template <class RegionalValues>
    void updateGpMaintTargetForGroups(const Group& group,
                                      const RegionalValues& regional_values,
                                      const double dt);

    /// update the number of wells that are actively under group control for a given group with name given by
    /// group_name its main usage is to detect cases where there is no wells under group control
    int updateGroupControlledWells(const bool is_production_group,
                                   const Phase injection_phase);

    void updateGroupProductionRates(const Group& group);

    void updateGroupTargetReduction(const Group& group,
                                    const bool is_injector);

    void updateNetworkLeafNodeProductionRates();

    void updateREINForGroups(const Group& group, bool sum_rank);

    void updateReservoirRatesInjectionGroups(const Group& group);

#ifdef RESERVOIR_COUPLING_ENABLED
    /// @brief Update the slave's GroupState cmodes from the master's active cmodes.
    ///
    /// For each slave group with a master-imposed target, sets the GroupState
    /// production/injection control mode to match the master's cmode. This ensures
    /// that all downstream consumers (constraint checks, guide rate fractions,
    /// well equation assembly) evaluate the correct rate type.
    ///
    /// The update is skipped when the GRUPSLAV filter flag is SLAV, meaning
    /// the slave ignores the master's control for that rate type.
    void updateSlaveGroupCmodesFromMaster();
#endif

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
                       const Group::ProductionCMode& offended_control) const;

private:
#ifdef RESERVOIR_COUPLING_ENABLED
    /// @brief Convert active phase index to ReservoirCoupling::Phase enum
    /// @param phase_pos Active phase index (0, 1, or 2 in a 3-phase model)
    /// @return The corresponding ReservoirCoupling::Phase enum value
    /// @note This uses the canonical phase ordering (Oil=0, Gas=1, Water=2)
    ReservoirCoupling::Phase activePhaseIdxToRescoupPhase_(int phase_pos) const;
#endif

    //! \brief Calculate group target by applying local rate adjustments and guide rate fractions through group hierarchy.
    //!
    //! This method encapsulates the core algorithm for calculating the portion of a top-level
    //! group target that should be assigned to a bottom group (well or sub-group). The algorithm:
    //! 1. Iterates through the group chain from top to bottom
    //! 2. At each level with a guide rate (or at the top level), subtracts local reductions (the amount of rate that is not available for the group to control)
    //! 3. At the local_reduction_level, adds back the bottom group's rate (if do_addback is true)
    //! 4. Multiplies by the guide rate fraction for each level
    //!
    //! \tparam ReductionLambda Callable type: (const std::string& group_name) -> Scalar
    //! \tparam FractionLambda Callable type: (const std::string& child) -> Scalar
    //! \param chain Group chain from control group (top) to bottom group
    //! \param orig_target Original target from top-level controlling group
    //! \param current_rate_available Current rate of the bottom group (for add-back)
    //! \param local_reduction_level Level where local reduction/add-back is applied
    //! \param is_production_group True for production, false for injection
    //! \param injection_phase Phase for injection groups (ignored for production)
    //! \param local_reduction_lambda Functor to get local reduction for a group
    //! \param local_fraction_lambda Functor to get local fraction for a child group
    //! \param do_addback Whether to perform add-back at local_reduction_level
    //! \return Target after applying reductions and fractions (before efficiency factor division)
    template<typename ReductionLambda, typename FractionLambda>
    Scalar applyReductionsAndFractions_(const std::vector<std::string>& chain,
                                        Scalar orig_target,
                                        Scalar current_rate_available,
                                        std::size_t local_reduction_level,
                                        bool is_production_group,
                                        Phase injection_phase,
                                        ReductionLambda&& local_reduction_lambda,
                                        FractionLambda&& local_fraction_lambda,
                                        bool do_addback) const;

    std::pair<Group::ProductionCMode, Scalar>
    checkProductionRateConstraint_(const Group& group,
                                   Group::ProductionCMode cmode,
                                   Group::ProductionCMode currentControl,
                                   Scalar target,
                                   Scalar current_rate) const;

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

    Scalar getInjectionGroupTargetForMode_(const Group& group,
        const Phase& injection_phase,
        const std::vector<Scalar>& resv_coeff,
        Group::InjectionCMode cmode) const;

    //! \brief Find the local reduction level in a group chain.
    //!
    //! The local reduction level is the deepest level in the chain (starting from level 1)
    //! where a group has both a guide rate and group-controlled wells (GCW > 0).
    //! This level determines where reductions are applied and add-back occurs.
    //!
    //! \param chain Group chain from control group (top) to bottom group
    //! \param is_production_group True for production, false for injection
    //! \param injection_phase Phase for injection groups (ignored for production)
    //! \return The local reduction level (0 if no intermediate group qualifies)
    std::size_t getLocalReductionLevel_(const std::vector<std::string>& chain,
        bool is_production_group,
        Phase injection_phase) const;

#ifdef RESERVOIR_COUPLING_ENABLED
    ReservoirCoupling::GrupSlav::FilterFlag getInjectionFilterFlag_(const std::string& group_name,
                                                                    Phase injection_phase) const;

    ReservoirCoupling::GrupSlav::FilterFlag getProductionFilterFlag_(
        const std::string& group_name,
        Group::ProductionCMode cmode) const;

    Scalar getReservoirCouplingMasterGroupRate_(const Group& group,
                                                const int phase_pos,
                                                ReservoirCoupling::RateKind kind) const;
#endif

    Scalar getProductionConstraintTarget_(const Group& group,
                                          Group::ProductionCMode cmode,
                                          const Group::ProductionControls& controls) const;

    Scalar getProductionGroupTargetForMode_(const Group& group, Group::ProductionCMode cmode) const;

    Scalar getSatelliteRate_(const Group& group,
        const int phase_pos,
        const bool res_rates,
        const bool is_injector) const;


    /// Check if a production auto choke group is underperforming its target rate.
    /// Returns true if the group's current rate is below its allocated target,
    /// which means wells should be excluded from the GCW count.
    bool isAutoChokeGroupUnderperforming_(const Group& group) const;

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

    Scalar sumProductionRate_(const Group& group, Group::ProductionCMode cmode) const;

    int updateGroupControlledWellsRecursive_(const std::string& group_name,
                                             const bool is_production_group,
                                             const Phase injection_phase);

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
    bool terminal_output_ {false};
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
    if (this->isReservoirCouplingMasterGroup(group)) {
        // GPMAINT is not supported for reservoir coupling master groups since master groups do not have
        //   subordinate wells in the master reservoir, so the slaves cannot influence the master reservoir's
        //   average pressure.
        //   Even specifying GPMAINT on a group superior to the master group might not make sense, since if the
        //   superior target is distributed down to the master group with guide rate fractions, adjusting
        //   the master group's target (that is sent to the slave) could only indirectly influence the master
        //   reservoir's average pressure by affecting the guide rate fractions distributed to actual wells
        //   in the master reservoir.
        OPM_DEFLOG_THROW(
            std::runtime_error,
            "GPMAINT is not supported for reservoir coupling master groups.",
            this->deferredLogger()
        );
        return;
    }
    else if (this->isReservoirCouplingSlaveGroup(group)) {
        // GPMAINT is not supported for reservoir coupling slave groups since their targets will be overridden
        //   by the corresponding master group's target anyway.
        OPM_DEFLOG_THROW(
            std::runtime_error,
            "GPMAINT is not supported for reservoir coupling slave groups.",
            this->deferredLogger()
        );
        return;
    }
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
