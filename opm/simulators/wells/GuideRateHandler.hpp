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
#include <opm/simulators/utils/BlackoilPhases.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/WellState.hpp>
namespace Opm {

template<class Scalar>
class GuideRateHandler {
public:
    using Potentials = ReservoirCoupling::Potentials;

    class UpdateGuideRates {
    public:
        UpdateGuideRates(
            GuideRateHandler<Scalar> &parent,
            const int report_step_idx,
            const double sim_time,
            const WellState<Scalar> &well_state,
            const GroupState<Scalar> &group_state,
            const int num_phases
        );

        const Parallel::Communication &comm() const { return this->parent_.comm_; }
        DeferredLogger &deferredLogger() { return this->parent_.deferredLogger(); }
        GuideRate &guideRate() { return this->parent_.guide_rate_; }
        bool isReservoirCouplingMaster() const { return this->parent_.isReservoirCouplingMaster(); }
        const PhaseUsage &phaseUsage() const { return this->parent_.phase_usage_; }
        ReservoirCouplingMaster& reservoirCouplingMaster() {
             return this->parent_.reservoirCouplingMaster();
        }
        const SummaryState &summaryState() const { return this->parent_.summary_state_; }
        const Schedule &schedule() const { return this->parent_.schedule_; }
        void update();

    private:
#ifdef RESERVOIR_COUPLING_ENABLED
        bool isMasterGroup_(const Group& group);
#endif
        void updateGuideRatesForInjectionGroups_(const Group& group);
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
        const GroupState<Scalar> &group_state_;
        const int num_phases_;
        const UnitSystem& unit_system_;
    };

    GuideRateHandler(
        const Schedule& schedule,
        const PhaseUsage& phase_usage,
        const SummaryState& summary_state,
        const Parallel::Communication& comm,
        GuideRate& guide_rate
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
    void sendSlaveGroupPotentialsToMaster();
    void setReservoirCouplingMaster(ReservoirCouplingMaster *reservoir_coupling_master) {
        this->reservoir_coupling_master_ = reservoir_coupling_master;
    }
    void setReservoirCouplingSlave(ReservoirCouplingSlave *reservoir_coupling_slave) {
        this->reservoir_coupling_slave_ = reservoir_coupling_slave;
    }
#endif
    DeferredLogger& deferredLogger();
    void setLogger(DeferredLogger *deferred_logger);
    void updateGuideRates(
        const int report_step_idx,
        const double sim_time,
        WellState<Scalar> &well_state,
        GroupState<Scalar> &group_state
    );
private:
    const Schedule& schedule_;
    const PhaseUsage& phase_usage_;
    const SummaryState& summary_state_;
    const Parallel::Communication& comm_;
    GuideRate& guide_rate_;
    DeferredLogger *deferred_logger_ = nullptr;
#ifdef RESERVOIR_COUPLING_ENABLED
    ReservoirCouplingMaster *reservoir_coupling_master_ = nullptr;
    ReservoirCouplingSlave *reservoir_coupling_slave_ = nullptr;
#endif
};

} // namespace Opm

#endif // OPM_GUIDERATE_HANDLER_HPP
