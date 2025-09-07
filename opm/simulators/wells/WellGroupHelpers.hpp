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

#include "opm/material/fluidsystems/PhaseUsageInfo.hpp"

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
class Schedule;
template<class Scalar> class VFPProdProperties;
template<typename Scalar, typename IndexTraits> class WellState;
class FieldPropsManager;

namespace Network { class ExtNetwork; }

template<typename Scalar, typename IndexTraits>
class WellGroupHelpers
{
public:

    using WellStateType = WellState<Scalar, IndexTraits>;

    static Scalar sumWellPhaseRates(bool res_rates,
                                    const Opm::Group& group,
                                    const Opm::Schedule& schedule,
                                    const WellStateType& wellState,
                                    const int reportStepIdx,
                                    const int phasePos,
                                    const bool injector,
                                    const bool network = false);

    static Scalar satelliteInjectionRate(const ScheduleState& sched,
                                         const Group& group,
                                         const PhaseUsageInfo<IndexTraits>& pu,
                                         const int phase_pos,
                                         bool res_rates);

    static Scalar satelliteProductionRate(const ScheduleState& sched,
                                          const Group& group, 
                                          const GSatProd::GSatProdGroup::Rate rateComp, 
                                          bool res_rates);

    static std::optional<GSatProd::GSatProdGroup::Rate>
    selectRateComponent(const PhaseUsageInfo<IndexTraits>& pu, const int phasePos);

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
                                      const WellStateType& wellState,
                                      const int reportStepIdx,
                                      const int phasePos,
                                      const bool injector);

    /// Returns the name of the worst offending well and its fraction (i.e. violated_phase / preferred_phase)
    static std::pair<std::optional<std::string>, Scalar>
    worstOffendingWell(const Group& group,
                       const Schedule& schedule,
                       const int reportStepIdx,
                       const Group::ProductionCMode& offendedControl,
                       const Parallel::Communication& comm,
                       const WellStateType& wellState,
                       DeferredLogger& deferred_logger);

    static Scalar sumWellResRates(const Group& group,
                                  const Schedule& schedule,
                                  const WellStateType& wellState,
                                  const int reportStepIdx,
                                  const int phasePos,
                                  const bool injector);

    static Scalar sumSolventRates(const Group& group,
                                  const Schedule& schedule,
                                  const WellStateType& wellState,
                                  const int reportStepIdx,
                                  const bool injector);

    static void updateGroupTargetReduction(const Group& group,
                                           const Schedule& schedule,
                                           const int reportStepIdx,
                                           const bool isInjector,
                                           const GuideRate& guide_rate,
                                           const WellStateType& wellState,
                                           const SummaryState& summaryState,
                                           GroupState<Scalar>& group_state,
                                           std::vector<Scalar>& groupTargetReduction);

    static void updateVREPForGroups(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const WellStateType& wellState,
                                    GroupState<Scalar>& group_state);

    template <class RegionalValues>
    static void updateGpMaintTargetForGroups(const Group& group,
                                             const Schedule& schedule,
                                             const RegionalValues& regional_values,
                                             const int reportStepIdx,
                                             const double dt,
                                             const WellStateType& well_state,
                                             GroupState<Scalar>& group_state);

    static void updateReservoirRatesInjectionGroups(const Group& group,
                                                    const Schedule& schedule,
                                                    const int reportStepIdx,
                                                    const WellStateType& wellState,
                                                    GroupState<Scalar>& group_state);

    static void updateSurfaceRatesInjectionGroups(const Group& group,
                                                  const Schedule& schedule,
                                                  const int reportStepIdx,
                                                  const WellStateType& wellState,
                                                  GroupState<Scalar>& group_state);

    static void updateWellRates(const Group& group,
                                const Schedule& schedule,
                                const int reportStepIdx,
                                const WellStateType& wellStateNupcol,
                                WellStateType& wellState);

    static void updateGroupProductionRates(const Group& group,
                                           const Schedule& schedule,
                                           const int reportStepIdx,
                                           const WellStateType& wellState,
                                           GroupState<Scalar>& group_state);

    static void updateNetworkLeafNodeProductionRates(const Schedule& schedule,
                                                     const int reportStepIdx,
                                                     const WellStateType& wellState,
                                                     GroupState<Scalar>& group_state);


    static void updateWellRatesFromGroupTargetScale(const Scalar scale,
                                                    const Group& group,
                                                    const Schedule& schedule,
                                                    const int reportStepIdx,
                                                    bool isInjector,
                                                    const GroupState<Scalar>& group_state,
                                                    WellStateType& wellState);

    static void updateREINForGroups(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const SummaryState& st,
                                    const WellStateType& wellState,
                                    GroupState<Scalar>& group_state,
                                    bool sum_rank);


    static std::map<std::string, Scalar>
    computeNetworkPressures(const Network::ExtNetwork& network,
                            const WellStateType& well_state,
                            const GroupState<Scalar>& group_state,
                            const VFPProdProperties<Scalar>& vfp_prod_props,
                            const Schedule& schedule,
                            const Parallel::Communication& comm,
                            const int report_time_step);

    static GuideRate::RateVector
    getWellRateVector(const WellStateType& well_state,
                      const std::string& name);

    static GuideRate::RateVector
    getProductionGroupRateVector(const GroupState<Scalar>& group_state,
                                 const PhaseUsageInfo<IndexTraits>& pu,
                                 const std::string& group_name);

    static Scalar getGuideRate(const std::string& name,
                               const Schedule& schedule,
                               const WellStateType& wellState,
                               const GroupState<Scalar>& group_state,
                               const int reportStepIdx,
                               const GuideRate* guideRate,
                               const GuideRateModel::Target target);

    static Scalar getGuideRateInj(const std::string& name,
                                  const Schedule& schedule,
                                  const WellStateType& wellState,
                                  const GroupState<Scalar>& group_state,
                                  const int reportStepIdx,
                                  const GuideRate* guideRate,
                                  const GuideRateModel::Target target,
                                  const Phase& injectionPhase);

    /// update the number of wells that are actively under group control for a given group with name given by group_name
    /// its main usage is to detect cases where there is no wells under group control
    static int updateGroupControlledWells(const Schedule& schedule,
                                          const WellStateType& well_state,
                                          GroupState<Scalar>& group_state,
                                          const SummaryState& summary_state,
                                          const GuideRate* guideRate,
                                          const int report_step,
                                          const std::string& group_name,
                                          const bool is_production_group,
                                          const Phase injection_phase);

    /// returns the number of wells that are actively under group control for a given group with name given by group_name
    static int groupControlledWells(const Schedule& schedule,
                                    const WellStateType& well_state,
                                    const GroupState<Scalar>& group_state,
                                    const int report_step,
                                    const std::string& group_name,
                                    const std::string& always_included_child,
                                    const bool is_production_group,
                                    const Phase injection_phase);

    static std::pair<bool, Scalar>
    checkGroupConstraintsInj(const std::string& name,
                             const std::string& parent,
                             const Group& group,
                             const WellStateType& wellState,
                             const GroupState<Scalar>& group_state,
                             const int reportStepIdx,
                             const GuideRate* guideRate,
                             const Scalar* rates,
                             Phase injectionPhase,
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
                               const WellStateType& wellState,
                               const GroupState<Scalar>& group_state,
                               const int reportStepIdx,
                               const GuideRate* guideRate,
                               const Scalar* rates,
                               Phase injectionPhase,
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

    
    // check if well/group bottom is a sub well/group of the group top
    static bool
    isInGroupChainTopBot(const std::string& bottom,
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
                              const WellStateType& wellState,
                              const GroupState<Scalar>& group_state,
                              const int reportStepIdx,
                              const GuideRate* guideRate,
                              const Scalar* rates,
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
                               const WellStateType& wellState,
                               const GroupState<Scalar>& group_state,
                               const int reportStepIdx,
                               const GuideRate* guideRate,
                               const Scalar* rates,
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
                                                   std::map<std::string, std::unique_ptr<AverageRegionalPressureType>>& regionalAveragePressureCalculator);
};

} // namespace Opm

#endif
