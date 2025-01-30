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

#ifndef OPM_GUIDERATE_HANDLER_HPP
#define OPM_GUIDERATE_HANDLER_HPP

#include <opm/simulators/wells/rescoup/RescoupProxy.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/output/data/GuideRateValue.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/BlackoilWellModelGuideRates.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/WellState.hpp>
namespace Opm {

/**
 * @brief Handles computation and reporting of guide rates for wells and groups.
 *
 * This class manages the update and dumping of guide rates, which are used to
 * control and monitor the production and injection targets for wells and groups.
 *
 * It integrates with the well model to compute guide rates and supports optional
 * reservoir coupling to exchange guide rate information between master and slave models.
 *
 * @tparam Scalar The scalar type (e.g., float or double) used in computations.
 */
template<typename Scalar, typename IndexTraits>
class GuideRateHandler {
public:

#ifdef RESERVOIR_COUPLING_ENABLED
    using Potentials = ReservoirCoupling::Potentials<Scalar>;
#endif

    /**
     * @brief Responsible for formatting and printing guide rate information to logs.
     *
     * Used primarily for debugging and human-readable output. Dumps well and group
     * guide rates at a given report step and simulation time.
     */
    class GuideRateDumper {
    public:
        GuideRateDumper(
            GuideRateHandler<Scalar, IndexTraits> &parent, const int report_step_idx, const double sim_time
        );

        DeferredLogger &deferredLogger() { return this->parent_.deferredLogger(); }
        /**
         * @brief Dumps guide rates for all wells and groups in a hierarchical structure.
         */
        void dumpGuideRates();
    private:
        void dumpGuideRatesRecursive_(const Group& group, int level);
        void getGroupGuideRatesInjection_(
            const Group& group,
            const data::GroupGuideRates& group_guide_rate,
            std::vector<std::string>& msg_items
        ) const;
        void getGroupGuideRatesProduction_(
            const Group& group,
            const data::GroupGuideRates& group_guide_rate,
            std::vector<std::string>& msg_items
        ) const;
        /**
         * @brief Helper to print formatted group guide rate values.
         *
         * @param group The group to print guide rate info for.
         * @param level Nesting level in the group hierarchy.
         */
        void printGroupGuideRates_(const Group& group, int level);
        void printHeader_();
        void printFooter_();
        void printWellGuideRates_(const Well& well, int level);

        GuideRateHandler<Scalar, IndexTraits> &parent_;
        const int report_step_idx_;
        const double sim_time_;
        const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model_;
        const Schedule& schedule_;
        const Parallel::Communication& comm_;
        std::unordered_map<std::string, data::GuideRateValue> well_guide_rates_;
        std::unordered_map<std::string, data::GroupGuideRates> group_guide_rates_;
    };

    /**
     * @brief Computes and updates guide rate values for wells and groups.
     *
     * Updates are done for a specific report step and simulation time. This class
     * also handles computation of group potentials and supports reservoir coupling
     * features if enabled.
     */
    class UpdateGuideRates {
    public:
        UpdateGuideRates(
            GuideRateHandler<Scalar, IndexTraits> &parent,
            const int report_step_idx,
            const double sim_time,
            const WellState<Scalar, IndexTraits> &well_state,
            GroupState<Scalar> &group_state,
            const int num_phases
        );

        bool isReservoirCouplingMaster() const { return this->parent_.isReservoirCouplingMaster(); }
        ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster() {
            return this->parent_.reservoirCouplingMaster();
        }
        const Parallel::Communication &comm() const { return this->parent_.comm_; }
        DeferredLogger &deferredLogger() { return this->parent_.deferredLogger(); }
        GuideRate &guideRate() { return this->parent_.guide_rate_; }
        const PhaseUsageInfo<IndexTraits>& phaseUsage() const { return this->parent_.wellModel().phaseUsage(); }
        const SummaryState &summaryState() const { return this->parent_.summary_state_; }
        const Schedule &schedule() const { return this->parent_.schedule_; }
        /**
         * @brief Triggers the guide rate update process for the current simulation step.
         */
        void update();

    private:
#ifdef RESERVOIR_COUPLING_ENABLED
        bool isMasterGroup_(const Group& group);
#endif
        void updateGuideRatesForInjectionGroups_(const Group& group);
        /**
         * @brief Updates guide rates for all production groups recursively.
         *
         * @param group The root group to update.
         * @param pot Output parameter for the computed production potentials.
         */
        void updateGuideRatesForProductionGroups_(const Group& group, std::vector<Scalar>& pot);
        void updateGuideRatesForWells_();
#ifdef RESERVOIR_COUPLING_ENABLED
        void updateProductionGroupPotentialFromSlaveGroup_(
            const Group& group, std::vector<Scalar>& pot);
#endif
        void updateProductionGroupPotentialFromSubGroups(
            const Group& group, std::vector<Scalar>& pot);

        GuideRateHandler<Scalar, IndexTraits> &parent_;
        const int report_step_idx_;
        const double sim_time_;
        const WellState<Scalar, IndexTraits> &well_state_;
        GroupState<Scalar> &group_state_;
        const int num_phases_;
        const UnitSystem& unit_system_;
    };

    GuideRateHandler(
        BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model,
        const Schedule& schedule,
        const SummaryState& summary_state,
        const Parallel::Communication& comm
    );

    // === Reservoir Coupling ===

    /// @brief Get the reservoir coupling proxy
    ReservoirCoupling::Proxy<Scalar>& rescoup() { return rescoup_; }
    const ReservoirCoupling::Proxy<Scalar>& rescoup() const { return rescoup_; }

    bool isReservoirCouplingMaster() const { return rescoup_.isMaster(); }
    bool isReservoirCouplingSlave() const { return rescoup_.isSlave(); }

    ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster() { return rescoup_.master(); }
    ReservoirCouplingSlave<Scalar>& reservoirCouplingSlave() { return rescoup_.slave(); }

#ifdef RESERVOIR_COUPLING_ENABLED
    void sendSlaveGroupPotentialsToMaster(const GroupState<Scalar>& group_state);
    void setReservoirCouplingMaster(ReservoirCouplingMaster<Scalar>* master) {
        rescoup_.setMaster(master);
    }
    void setReservoirCouplingSlave(ReservoirCouplingSlave<Scalar>* slave) {
        rescoup_.setSlave(slave);
    }
#endif
    DeferredLogger& deferredLogger();
    /**
     * @brief Dumps guide rate information to the logger in a readable format.
     *
     * Used mainly for debugging or inspection during simulation development.
     *
     * @param report_step_idx Index of the current report step.
     * @param sim_time Current simulation time.
     */
    void debugDumpGuideRates(const int report_step_idx, const double sim_time);
    const Parallel::Communication& getComm() const { return comm_; }
    void setLogger(DeferredLogger *deferred_logger) { deferred_logger_ = deferred_logger; }
    const GuideRate& guideRate() { return guide_rate_; }
    const PhaseUsageInfo<IndexTraits>& phaseUsage() const { return well_model_.phaseUsage(); }
    const Schedule& schedule() const { return schedule_; }
    const SummaryState& summaryState() const { return summary_state_; }
    /**
     * @brief Updates guide rates for the current simulation step.
     *
     * @param report_step_idx Index of the current report step.
     * @param sim_time Simulation time at the current step.
     * @param well_state Well state object containing current well potentials.
     * @param group_state Group state object to update with computed guide rates.
     */
    void updateGuideRates(const int report_step_idx,
                          const double sim_time,
                          const WellState<Scalar, IndexTraits>& well_state,
                          GroupState<Scalar>& group_state);

    const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel() const { return well_model_; }
    BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel() { return well_model_; }
private:
    void debugDumpGuideRatesRecursive_(const Group& group) const;
    BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    const Parallel::Communication& comm_;
    GuideRate& guide_rate_;
    DeferredLogger *deferred_logger_ = nullptr;
    ReservoirCoupling::Proxy<Scalar> rescoup_{};
};

} // namespace Opm

#endif // OPM_GUIDERATE_HANDLER_HPP
