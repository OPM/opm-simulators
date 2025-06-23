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

#if HAVE_MPI
#define RESERVOIR_COUPLING_ENABLED
#endif
#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>
#endif
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/output/data/GuideRateValue.hpp>
#include <opm/simulators/utils/BlackoilPhases.hpp>
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
template<class Scalar>
class GuideRateHandler {
public:
#ifdef RESERVOIR_COUPLING_ENABLED
    using Potentials = ReservoirCoupling::Potentials;
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
            GuideRateHandler<Scalar> &parent, const int report_step_idx, const double sim_time
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

        GuideRateHandler<Scalar> &parent_;
        const int report_step_idx_;
        const double sim_time_;
        const BlackoilWellModelGeneric<Scalar>& well_model_;
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
            GuideRateHandler<Scalar> &parent,
            const int report_step_idx,
            const double sim_time,
            const WellState<Scalar> &well_state,
            GroupState<Scalar> &group_state,
            const int num_phases
        );

#ifdef RESERVOIR_COUPLING_ENABLED
        bool isReservoirCouplingMaster() const { return this->parent_.isReservoirCouplingMaster(); }
        ReservoirCouplingMaster& reservoirCouplingMaster() {
            return this->parent_.reservoirCouplingMaster();
        }
#endif
       const Parallel::Communication &comm() const { return this->parent_.comm_; }
        DeferredLogger &deferredLogger() { return this->parent_.deferredLogger(); }
        GuideRate &guideRate() { return this->parent_.guide_rate_; }
        const PhaseUsage &phaseUsage() const { return this->parent_.phase_usage_; }
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

        GuideRateHandler<Scalar> &parent_;
        const int report_step_idx_;
        const double sim_time_;
        const WellState<Scalar> &well_state_;
        GroupState<Scalar> &group_state_;
        const int num_phases_;
        const UnitSystem& unit_system_;
    };

    GuideRateHandler(
        BlackoilWellModelGeneric<Scalar>& well_model,
        const Schedule& schedule,
        const SummaryState& summary_state,
        const Parallel::Communication& comm
    );

#ifdef RESERVOIR_COUPLING_ENABLED
    bool isReservoirCouplingMaster() const {
        return this->reservoir_coupling_master_ != nullptr;
    }
    bool isReservoirCouplingSlave() const {
        return this->reservoir_coupling_slave_ != nullptr;
    }
    void receiveMasterGroupPotentialsFromSlaves();
    ReservoirCouplingMaster& reservoirCouplingMaster() { return *(this->reservoir_coupling_master_); }
    ReservoirCouplingSlave& reservoirCouplingSlave() { return *(this->reservoir_coupling_slave_); }
    void sendSlaveGroupPotentialsToMaster(const GroupState<Scalar>& group_state);
    void setReservoirCouplingMaster(ReservoirCouplingMaster *reservoir_coupling_master) {
        this->reservoir_coupling_master_ = reservoir_coupling_master;
    }
    void setReservoirCouplingSlave(ReservoirCouplingSlave *reservoir_coupling_slave) {
        this->reservoir_coupling_slave_ = reservoir_coupling_slave;
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
    void setLogger(DeferredLogger *deferred_logger);
    const Schedule& schedule() const { return schedule_; }
    /**
     * @brief Updates guide rates for the current simulation step.
     *
     * @param report_step_idx Index of the current report step.
     * @param sim_time Simulation time at the current step.
     * @param well_state Well state object containing current well potentials.
     * @param group_state Group state object to update with computed guide rates.
     */
    void updateGuideRates(
        const int report_step_idx,
        const double sim_time,
        WellState<Scalar> &well_state,
        GroupState<Scalar> &group_state
    );
    const BlackoilWellModelGeneric<Scalar>& wellModel() const { return well_model_; }
private:
    void debugDumpGuideRatesRecursive_(const Group& group) const;
    BlackoilWellModelGeneric<Scalar>& well_model_;
    const Schedule& schedule_;
    const SummaryState& summary_state_;
    const Parallel::Communication& comm_;
    const PhaseUsage& phase_usage_;
    GuideRate& guide_rate_;
    DeferredLogger *deferred_logger_ = nullptr;
#ifdef RESERVOIR_COUPLING_ENABLED
    ReservoirCouplingMaster *reservoir_coupling_master_ = nullptr;
    ReservoirCouplingSlave *reservoir_coupling_slave_ = nullptr;
#endif
};

} // namespace Opm

#endif // OPM_GUIDERATE_HANDLER_HPP
