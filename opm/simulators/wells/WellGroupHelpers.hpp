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

#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm {

class DeferredLogger;
class Group;
template<class Scalar> class GroupState;
namespace Network { class ExtNetwork; }
struct PhaseUsage;
class Schedule;
template<class Scalar> class VFPProdProperties;
template<class Scalar> class WellState;
class FieldPropsManager;

namespace Network { class ExtNetwork; }

template<class Scalar>
class WellGroupHelpers
{
public:

    static Scalar sumWellPhaseRates(bool res_rates,
                                    const Opm::Group& group,
                                    const Opm::Schedule& schedule,
                                    const Opm::WellState<Scalar>& wellState,
                                    const int reportStepIdx,
                                    const int phasePos,
                                    const bool injector,
                                    const bool network = false);

    static Scalar satelliteProduction(const ScheduleState& sched,
                                      const std::vector<std::string>& groups,
                                      const GSatProd::GSatProdGroup::Rate rateComp);

    static std::optional<GSatProd::GSatProdGroup::Rate>
    selectRateComponent(const PhaseUsage& pu, const int phasePos);

    static void setCmodeGroup(const Group& group,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              const int reportStepIdx,
                              GroupState<Scalar>& group_state);

    static void accumulateGroupEfficiencyFactor(const Group& group,
                                                const Schedule& schedule,
                                                const int reportStepIdx,
                                                Scalar& factor);

    static Scalar sumWellSurfaceRates(const Group& group,
                                      const Schedule& schedule,
                                      const WellState<Scalar>& wellState,
                                      const int reportStepIdx,
                                      const int phasePos,
                                      const bool injector);

    /// Returns the name of the worst offending well and its fraction (i.e. violated_phase / preferred_phase)
    static std::pair<std::optional<std::string>, Scalar>
    worstOffendingWell(const Group& group,
                       const Schedule& schedule,
                       const int reportStepIdx,
                       const Group::ProductionCMode& offendedControl,
                       const PhaseUsage& pu,
                       const Parallel::Communication& comm,
                       const WellState<Scalar>& wellState,
                       DeferredLogger& deferred_logger);

    static Scalar sumWellResRates(const Group& group,
                                  const Schedule& schedule,
                                  const WellState<Scalar>& wellState,
                                  const int reportStepIdx,
                                  const int phasePos,
                                  const bool injector);

    static Scalar sumSolventRates(const Group& group,
                                  const Schedule& schedule,
                                  const WellState<Scalar>& wellState,
                                  const int reportStepIdx,
                                  const bool injector);

    static void updateGroupTargetReduction(const Group& group,
                                           const Schedule& schedule,
                                           const int reportStepIdx,
                                           const bool isInjector,
                                           const PhaseUsage& pu,
                                           const GuideRate& guide_rate,
                                           const WellState<Scalar>& wellState,
                                           const SummaryState& summaryState,
                                           GroupState<Scalar>& group_state,
                                           std::vector<Scalar>& groupTargetReduction);

    static void updateGuideRates(const Group& group,
                                 const Schedule& schedule,
                                 const SummaryState& summary_state,
                                 const PhaseUsage& pu,
                                 int report_step,
                                 double sim_time,
                                 WellState<Scalar>& well_state,
                                 const GroupState<Scalar>& group_state,
                                 const Parallel::Communication& comm,
                                 GuideRate* guide_rate,
                                 std::vector<Scalar>& pot,
                                 DeferredLogger& deferred_logger);

    static void updateGuideRateForProductionGroups(const Group& group,
                                                   const Schedule& schedule,
                                                   const PhaseUsage& pu,
                                                   const int reportStepIdx,
                                                   const double& simTime,
                                                   WellState<Scalar>& wellState,
                                                   const GroupState<Scalar>& group_state,
                                                   const Parallel::Communication& comm,
                                                   GuideRate* guideRate,
                                                   std::vector<Scalar>& pot);

    static void updateGuideRatesForWells(const Schedule& schedule,
                                         const PhaseUsage& pu,
                                         const int reportStepIdx,
                                         const double& simTime,
                                         const WellState<Scalar>& wellState,
                                         const Parallel::Communication& comm,
                                         GuideRate* guideRate);

    static void updateGuideRatesForInjectionGroups(const Group& group,
                                                   const Schedule& schedule,
                                                   const SummaryState& summaryState,
                                                   const PhaseUsage& pu,
                                                   const int reportStepIdx,
                                                   const WellState<Scalar>& wellState,
                                                   const GroupState<Scalar>& group_state,
                                                   GuideRate* guideRate,
                                                   DeferredLogger& deferred_logger);

    static void updateVREPForGroups(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const WellState<Scalar>& wellState,
                                    GroupState<Scalar>& group_state);

    template <class RegionalValues>
    static void updateGpMaintTargetForGroups(const Group& group,
                                             const Schedule& schedule,
                                             const RegionalValues& regional_values,
                                             const int reportStepIdx,
                                             const double dt,
                                             const WellState<Scalar>& well_state,
                                             GroupState<Scalar>& group_state);

    static void updateReservoirRatesInjectionGroups(const Group& group,
                                                    const Schedule& schedule,
                                                    const int reportStepIdx,
                                                    const WellState<Scalar>& wellState,
                                                    GroupState<Scalar>& group_state);

    static void updateSurfaceRatesInjectionGroups(const Group& group,
                                                  const Schedule& schedule,
                                                  const int reportStepIdx,
                                                  const WellState<Scalar>& wellState,
                                                  GroupState<Scalar>& group_state);

    static void updateWellRates(const Group& group,
                                const Schedule& schedule,
                                const int reportStepIdx,
                                const WellState<Scalar>& wellStateNupcol,
                                WellState<Scalar>& wellState);

    static void updateGroupProductionRates(const Group& group,
                                           const Schedule& schedule,
                                           const int reportStepIdx,
                                           const WellState<Scalar>& wellState,
                                           GroupState<Scalar>& group_state);

    static void updateNetworkLeafNodeProductionRates(const Schedule& schedule,
                                                     const int reportStepIdx,
                                                     const WellState<Scalar>& wellState,
                                                     GroupState<Scalar>& group_state);


    static void updateWellRatesFromGroupTargetScale(const Scalar scale,
                                                    const Group& group,
                                                    const Schedule& schedule,
                                                    const int reportStepIdx,
                                                    bool isInjector,
                                                    const GroupState<Scalar>& group_state,
                                                    WellState<Scalar>& wellState);

    static void updateREINForGroups(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const PhaseUsage& pu,
                                    const SummaryState& st,
                                    const WellState<Scalar>& wellState,
                                    GroupState<Scalar>& group_state,
                                    bool sum_rank);


    static std::map<std::string, Scalar>
    computeNetworkPressures(const Network::ExtNetwork& network,
                            const WellState<Scalar>& well_state,
                            const GroupState<Scalar>& group_state,
                            const VFPProdProperties<Scalar>& vfp_prod_props,
                            const Schedule& schedule,
                            const Parallel::Communication& comm,
                            const int report_time_step);

    static GuideRate::RateVector
    getWellRateVector(const WellState<Scalar>& well_state,
                      const PhaseUsage& pu,
                      const std::string& name);

    static GuideRate::RateVector
    getProductionGroupRateVector(const GroupState<Scalar>& group_state,
                                 const PhaseUsage& pu,
                                 const std::string& group_name);

    static Scalar getGuideRate(const std::string& name,
                               const Schedule& schedule,
                               const WellState<Scalar>& wellState,
                               const GroupState<Scalar>& group_state,
                               const int reportStepIdx,
                               const GuideRate* guideRate,
                               const GuideRateModel::Target target,
                               const PhaseUsage& pu);

    static Scalar getGuideRateInj(const std::string& name,
                                  const Schedule& schedule,
                                  const WellState<Scalar>& wellState,
                                  const GroupState<Scalar>& group_state,
                                  const int reportStepIdx,
                                  const GuideRate* guideRate,
                                  const GuideRateModel::Target target,
                                  const Phase& injectionPhase,
                                  const PhaseUsage& pu);

    static int groupControlledWells(const Schedule& schedule,
                                    const WellState<Scalar>& well_state,
                                    const GroupState<Scalar>& group_state,
                                    const SummaryState& summary_state,
                                    const GuideRate* guideRate,
                                    const int report_step,
                                    const std::string& group_name,
                                    const std::string& always_included_child,
                                    const bool is_production_group,
                                    const Phase injection_phase);

    static std::pair<bool, Scalar>
    checkGroupConstraintsInj(const std::string& name,
                             const std::string& parent,
                             const Group& group,
                             const WellState<Scalar>& wellState,
                             const GroupState<Scalar>& group_state,
                             const int reportStepIdx,
                             const GuideRate* guideRate,
                             const Scalar* rates,
                             Phase injectionPhase,
                             const PhaseUsage& pu,
                             const Scalar efficiencyFactor,
                             const Schedule& schedule,
                             const SummaryState& summaryState,
                             const std::vector<Scalar>& resv_coeff,
                             const bool check_guide_rate,
                             DeferredLogger& deferred_logger);

    static Scalar
    getWellGroupTargetInjector(const std::string& name,
                               const std::string& parent,
                               const Group& group,
                               const WellState<Scalar>& wellState,
                               const GroupState<Scalar>& group_state,
                               const int reportStepIdx,
                               const GuideRate* guideRate,
                               const Scalar* rates,
                               Phase injectionPhase,
                               const PhaseUsage& pu,
                               const Scalar efficiencyFactor,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const std::vector<Scalar>& resv_coeff,
                               DeferredLogger& deferred_logger);

    static std::vector<std::string>
    groupChainTopBot(const std::string& bottom,
                     const std::string& top,
                     const Schedule& schedule,
                     const int report_step);

    static std::string
    control_group(const Group& group,
                  const GroupState<Scalar>& group_state,
                  const int reportStepIdx,
                  const Schedule& schedule);

    static std::pair<bool, Scalar>
    checkGroupConstraintsProd(const std::string& name,
                              const std::string& parent,
                              const Group& group,
                              const WellState<Scalar>& wellState,
                              const GroupState<Scalar>& group_state,
                              const int reportStepIdx,
                              const GuideRate* guideRate,
                              const Scalar* rates,
                              const PhaseUsage& pu,
                              const Scalar efficiencyFactor,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              const std::vector<Scalar>& resv_coeff,
                              const bool check_guide_rate,
                              DeferredLogger& deferred_logger);
    static Scalar
    getWellGroupTargetProducer(const std::string& name,
                               const std::string& parent,
                               const Group& group,
                               const WellState<Scalar>& wellState,
                               const GroupState<Scalar>& group_state,
                               const int reportStepIdx,
                               const GuideRate* guideRate,
                               const Scalar* rates,
                               const PhaseUsage& pu,
                               const Scalar efficiencyFactor,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const std::vector<Scalar>& resv_coeff,
                               DeferredLogger& deferred_logger);

    template <class AverageRegionalPressureType>
    static void setRegionAveragePressureCalculator(const Group& group,
                                                   const Schedule& schedule,
                                                   const int reportStepIdx,
                                                   const FieldPropsManager& fp,
                                                   const PhaseUsage& pu,
                                                   std::map<std::string, std::unique_ptr<AverageRegionalPressureType>>& regionalAveragePressureCalculator);
};

} // namespace Opm

#endif
