/*
  Copyright 2019 Norce.

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


#ifndef OPM_WELLGROUPHELPERS_HEADER_INCLUDED
#define OPM_WELLGROUPHELPERS_HEADER_INCLUDED

#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm
{

class DeferredLogger;
class Group;
class GroupState;
namespace Network { class ExtNetwork; }
struct PhaseUsage;
class Schedule;
class VFPProdProperties;
class WellState;

template <typename>
class WellContainer;

namespace Network { class ExtNetwork; }

namespace WellGroupHelpers
{



    void setCmodeGroup(const Group& group,
                       const Schedule& schedule,
                       const SummaryState& summaryState,
                       const int reportStepIdx,
                       WellState& wellState,
                       GroupState& group_state);

    void accumulateGroupEfficiencyFactor(const Group& group,
                                         const Schedule& schedule,
                                         const int reportStepIdx,
                                         double& factor);

    double sumWellSurfaceRates(const Group& group,
                               const Schedule& schedule,
                               const WellState& wellState,
                               const int reportStepIdx,
                               const int phasePos,
                               const bool injector);

    double sumWellResRates(const Group& group,
                           const Schedule& schedule,
                           const WellState& wellState,
                           const int reportStepIdx,
                           const int phasePos,
                           const bool injector);

    double sumSolventRates(const Group& group,
                           const Schedule& schedule,
                           const WellState& wellState,
                           const int reportStepIdx,
                           const bool injector);

    void updateGroupTargetReduction(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const bool isInjector,
                                    const PhaseUsage& pu,
                                    const GuideRate& guide_rate,
                                    const WellState& wellStateNupcol,
                                    WellState& wellState,
                                    GroupState& group_state,
                                    std::vector<double>& groupTargetReduction);

    template <class Comm>
    void updateGuideRates(const Group& group,
                          const Schedule& schedule,
                          const SummaryState& summary_state,
                          const PhaseUsage& pu,
                          int report_step,
                          double sim_time,
                          WellState& well_state,
                          const GroupState& group_state,
                          const Comm& comm,
                          GuideRate* guide_rate,
                          std::vector<double>& pot,
                          Opm::DeferredLogger& deferred_logge);

    template <class Comm>
    void updateGuideRateForProductionGroups(const Group& group,
                                            const Schedule& schedule,
                                            const PhaseUsage& pu,
                                            const int reportStepIdx,
                                            const double& simTime,
                                            WellState& wellState,
                                            const GroupState& group_state,
                                            const Comm& comm,
                                            GuideRate* guideRate,
                                            std::vector<double>& pot);

    template <class Comm>
    void updateGuideRatesForWells(const Schedule& schedule,
                                  const PhaseUsage& pu,
                                  const int reportStepIdx,
                                  const double& simTime,
                                  const WellState& wellState,
                                  const Comm& comm,
                                  GuideRate* guideRate);

    void updateGuideRatesForInjectionGroups(const Group& group,
                                            const Schedule& schedule,
                                            const SummaryState& summaryState,
                                            const Opm::PhaseUsage& pu,
                                            const int reportStepIdx,
                                            const WellState& wellState,
                                            const GroupState& group_state,
                                            GuideRate* guideRate,
                                            Opm::DeferredLogger& deferred_logger);

    void updateVREPForGroups(const Group& group,
                             const Schedule& schedule,
                             const int reportStepIdx,
                             const WellState& wellStateNupcol,
                             WellState& wellState,
                             GroupState& group_state);

    void updateReservoirRatesInjectionGroups(const Group& group,
                                             const Schedule& schedule,
                                             const int reportStepIdx,
                                             const WellState& wellStateNupcol,
                                             WellState& wellState,
                                             GroupState& group_state);

    void updateWellRates(const Group& group,
                         const Schedule& schedule,
                         const int reportStepIdx,
                         const WellState& wellStateNupcol,
                         WellState& wellState);

    void updateGroupProductionRates(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const WellState& wellStateNupcol,
                                    WellState& wellState,
                                    GroupState& group_state);

    void updateWellRatesFromGroupTargetScale(const double scale,
                                             const Group& group,
                                             const Schedule& schedule,
                                             const int reportStepIdx,
                                             bool isInjector,
                                             const GroupState& group_state,
                                             WellState& wellState);

    void updateREINForGroups(const Group& group,
                             const Schedule& schedule,
                             const int reportStepIdx,
                             const PhaseUsage& pu,
                             const SummaryState& st,
                             const WellState& wellStateNupcol,
                             WellState& wellState,
                             GroupState& group_state);

    std::map<std::string, double>
    computeNetworkPressures(const Opm::Network::ExtNetwork& network,
                            const WellState& well_state,
                            const GroupState& group_state,
                            const VFPProdProperties& vfp_prod_props,
                            const Schedule& schedule,
                            const int report_time_step);

    GuideRate::RateVector
    getWellRateVector(const WellState& well_state, const PhaseUsage& pu, const std::string& name);

    GuideRate::RateVector
    getProductionGroupRateVector(const GroupState& group_state, const PhaseUsage& pu, const std::string& group_name);

    double getGuideRate(const std::string& name,
                        const Schedule& schedule,
                        const WellState& wellState,
                        const GroupState& group_state,
                        const int reportStepIdx,
                        const GuideRate* guideRate,
                        const GuideRateModel::Target target,
                        const PhaseUsage& pu);


    double getGuideRateInj(const std::string& name,
                           const Schedule& schedule,
                           const WellState& wellState,
                           const GroupState& group_state,
                           const int reportStepIdx,
                           const GuideRate* guideRate,
                           const GuideRateModel::Target target,
                           const Phase& injectionPhase,
                           const PhaseUsage& pu);

    int groupControlledWells(const Schedule& schedule,
                             const WellState& well_state,
                             const GroupState& group_state,
                             const int report_step,
                             const std::string& group_name,
                             const std::string& always_included_child,
                             const bool is_production_group,
                             const Phase injection_phase);


    class FractionCalculator
    {
    public:
        FractionCalculator(const Schedule& schedule,
                           const WellState& well_state,
                           const GroupState& group_state,
                           const int report_step,
                           const GuideRate* guide_rate,
                           const GuideRateModel::Target target,
                           const PhaseUsage& pu,
                           const bool is_producer,
                           const Phase injection_phase);
        double fraction(const std::string& name, const std::string& control_group_name, const bool always_include_this);
        double localFraction(const std::string& name, const std::string& always_included_child);

    private:
        std::string parent(const std::string& name);
        double guideRateSum(const Group& group, const std::string& always_included_child);
        double guideRate(const std::string& name, const std::string& always_included_child);
        int groupControlledWells(const std::string& group_name, const std::string& always_included_child);
        GuideRate::RateVector getGroupRateVector(const std::string& group_name);
        const Schedule& schedule_;
        const WellState& well_state_;
        const GroupState& group_state_;
        int report_step_;
        const GuideRate* guide_rate_;
        GuideRateModel::Target target_;
        const PhaseUsage& pu_;
        bool is_producer_;
        Phase injection_phase_;
    };


    std::pair<bool, double> checkGroupConstraintsInj(const std::string& name,
                                                     const std::string& parent,
                                                     const Group& group,
                                                     const WellState& wellState,
                                                     const GroupState& group_state,
                                                     const int reportStepIdx,
                                                     const GuideRate* guideRate,
                                                     const double* rates,
                                                     Phase injectionPhase,
                                                     const PhaseUsage& pu,
                                                     const double efficiencyFactor,
                                                     const Schedule& schedule,
                                                     const SummaryState& summaryState,
                                                     const std::vector<double>& resv_coeff,
                                                     DeferredLogger& deferred_logger);






    std::vector<std::string> groupChainTopBot(const std::string& bottom,
                                              const std::string& top,
                                              const Schedule& schedule,
                                              const int report_step);




    std::pair<bool, double> checkGroupConstraintsProd(const std::string& name,
                                                      const std::string& parent,
                                                      const Group& group,
                                                      const WellState& wellState,
                                                      const GroupState& group_state,
                                                      const int reportStepIdx,
                                                      const GuideRate* guideRate,
                                                      const double* rates,
                                                      const PhaseUsage& pu,
                                                      const double efficiencyFactor,
                                                      const Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      const std::vector<double>& resv_coeff,
                                                      DeferredLogger& deferred_logger);


} // namespace WellGroupHelpers

} // namespace Opm

#endif
