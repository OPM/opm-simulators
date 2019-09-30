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

#include <vector>

namespace Opm {


    namespace wellGroupHelpers
    {

    inline void setCmodeGroup(const Group2& group, const Schedule& schedule, const SummaryState& summaryState, const int reportStepIdx, WellStateFullyImplicitBlackoil& wellState) {

        for (const std::string& groupName : group.groups()) {
            setCmodeGroup( schedule.getGroup2(groupName, reportStepIdx), schedule, summaryState, reportStepIdx, wellState);
        }

        // use FLD as default control
        wellState.setCurrentInjectionGroupControl(group.name(), Group2::InjectionCMode::FLD);
        wellState.setCurrentProductionGroupControl(group.name(), Group2::ProductionCMode::FLD);

        if (group.isInjectionGroup()) {
            const auto controls = group.injectionControls(summaryState);
            wellState.setCurrentInjectionGroupControl(group.name(), controls.cmode);
        }
        if (group.isProductionGroup()) {
            const auto controls = group.productionControls(summaryState);
            wellState.setCurrentProductionGroupControl(group.name(), controls.cmode);
        }
    }


    inline void accumulateGroupEfficiencyFactor(const Group2& group, const Schedule& schedule, const int reportStepIdx, double& factor) {
        factor *= group.getGroupEfficiencyFactor();
        if (group.parent() != "FIELD")
            accumulateGroupEfficiencyFactor(schedule.getGroup2(group.parent(), reportStepIdx), schedule, reportStepIdx, factor);
    }

    inline void setGroupControl(const Group2& group, const Schedule& schedule, const int reportStepIdx, const bool injector, WellStateFullyImplicitBlackoil& wellState) {

        for (const std::string& groupName : group.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);
            setGroupControl(groupTmp, schedule, reportStepIdx, injector, wellState);
            if (injector)
                wellState.setCurrentInjectionGroupControl(groupName, Group2::InjectionCMode::FLD);
            else
                wellState.setCurrentProductionGroupControl(groupName, Group2::ProductionCMode::FLD);
        }

        const auto& end = wellState.wellMap().end();
        for (const std::string& wellName : group.wells()) {
            const auto& it = wellState.wellMap().find( wellName );
            if (it == end)  // the well is not found
                continue;

            int well_index = it->second[0];
            const auto& wellEcl = schedule.getWell2(wellName, reportStepIdx);

            if (wellEcl.getStatus() == Well2::Status::SHUT)
                continue;

            if (!wellEcl.isAvailableForGroupControl())
                continue;

            if (wellEcl.isProducer() && !injector)
                wellState.currentProductionControls()[well_index] = Well2::ProducerCMode::GRUP;

            if (wellEcl.isInjector() && injector)
                wellState.currentInjectionControls()[well_index] = Well2::InjectorCMode::GRUP;
        }
    }


    inline void computeGroupTargetReduction(const Group2& group, const WellStateFullyImplicitBlackoil& wellState, const Schedule& schedule, const int reportStepIdx, const int phasePos, const bool isInjector, double& groupTargetReduction )
    {

        for (const std::string& groupName : group.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);
            computeGroupTargetReduction(groupTmp, wellState, schedule, reportStepIdx, phasePos, isInjector, groupTargetReduction);
        }
        for (const std::string& wellName : group.wells()) {
            const auto& wellTmp = schedule.getWell2(wellName, reportStepIdx);

            if (wellTmp.isProducer() && isInjector)
                 continue;

            if (wellTmp.isInjector() && !isInjector)
                 continue;

            if (wellTmp.getStatus() == Well2::Status::SHUT)
                continue;

            const auto& end = wellState.wellMap().end();
                const auto& it = wellState.wellMap().find( wellName );
                if (it == end)  // the well is not found
                    continue;

            int well_index = it->second[0];
            const auto wellrate_index = well_index * wellState.numPhases();
            // add contributino from wells not under group control
            if (isInjector) {
                if (wellState.currentInjectionControls()[well_index] != Well2::InjectorCMode::GRUP)
                    groupTargetReduction += wellState.wellRates()[wellrate_index + phasePos];
            } else {
                if (wellState.currentProductionControls()[well_index] !=  Well2::ProducerCMode::GRUP)
                    groupTargetReduction -= wellState.wellRates()[wellrate_index + phasePos];
            }
        }
    }


    inline double sumWellPhaseRates(const std::vector<double>& rates, const Group2& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos,
                                    const bool injector) {

        double rate = 0.0;
        for (const std::string& groupName : group.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);
            rate += groupTmp.getGroupEfficiencyFactor()*sumWellPhaseRates(rates, groupTmp, schedule, wellState, reportStepIdx, phasePos, injector);
        }
        const auto& end = wellState.wellMap().end();
        for (const std::string& wellName : group.wells()) {
            const auto& it = wellState.wellMap().find( wellName );
            if (it == end)  // the well is not found
                continue;

            int well_index = it->second[0];

            const auto& wellEcl = schedule.getWell2(wellName, reportStepIdx);
            //only count producers or injectors
            if ( (wellEcl.isProducer() && injector) ||  (wellEcl.isInjector() && !injector))
                continue;

            if (wellEcl.getStatus() == Well2::Status::SHUT)
                continue;

            double factor = wellEcl.getEfficiencyFactor();
            const auto wellrate_index = well_index * wellState.numPhases();
            if (injector)
                rate += factor * rates[ wellrate_index + phasePos];
            else
                rate -= factor * rates[ wellrate_index + phasePos];
        }
        return rate;
    }

    inline double sumWellRates(const Group2& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos, const bool injector) {
        return sumWellPhaseRates(wellState.wellRates(), group, schedule, wellState, reportStepIdx, phasePos, injector);
    }

    inline double sumWellResRates(const Group2& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos, const bool injector) {
        return sumWellPhaseRates(wellState.wellReservoirRates(), group, schedule, wellState, reportStepIdx, phasePos, injector);
    }


    inline void updateGuideRateForGroups(const Group2& group, const Schedule& schedule, const PhaseUsage& pu, const int reportStepIdx, const double& simTime, GuideRate* guideRate, WellStateFullyImplicitBlackoil& wellState) {
        for (const std::string& groupName : group.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);
            updateGuideRateForGroups(groupTmp, schedule, pu, reportStepIdx, simTime, guideRate, wellState);
        }
        bool isInjector = group.isInjectionGroup();
        bool isProducer = group.isProductionGroup();

        if (!isInjector && !isProducer)
            return;

        double oilPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Liquid])
            oilPot = sumWellPhaseRates(wellState.wellPotentials(), group, schedule, wellState, reportStepIdx, pu.phase_pos[BlackoilPhases::Liquid], isInjector);

        double gasPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Vapour])
            gasPot = sumWellPhaseRates(wellState.wellPotentials(), group, schedule, wellState, reportStepIdx, pu.phase_pos[BlackoilPhases::Vapour], isInjector);

        double waterPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Aqua])
            waterPot = sumWellPhaseRates(wellState.wellPotentials(), group, schedule, wellState, reportStepIdx, pu.phase_pos[BlackoilPhases::Aqua], isInjector);

        guideRate->compute(group.name(), reportStepIdx, simTime, oilPot, gasPot, waterPot);
    }

    inline double wellFractionFromGuideRates(const Well2& well, const Schedule& schedule, const int reportStepIdx, const GuideRate* guideRate, const Well2::GuideRateTarget& wellTarget, const bool isInjector) {
        double groupTotalGuideRate = 0.0;
        const Group2& groupTmp = schedule.getGroup2(well.groupName(), reportStepIdx);
        for (const std::string& wellName : groupTmp.wells()) {
            const auto& wellTmp = schedule.getWell2(wellName, reportStepIdx);

            if (wellTmp.isProducer() && isInjector)
                 continue;

            if (wellTmp.isInjector() && !isInjector)
                 continue;

            if (wellTmp.getStatus() == Well2::Status::SHUT)
                continue;

            groupTotalGuideRate += guideRate->get(wellName, wellTarget);
        }
        double wellGuideRate = guideRate->get(well.name(), wellTarget);
        return wellGuideRate / groupTotalGuideRate;
    }

    inline double groupFractionFromGuideRates(const Group2& group, const Schedule& schedule, const int reportStepIdx, const GuideRate* guideRate, const Group2::GuideRateTarget& groupTarget, const bool isInjector) {
        double groupTotalGuideRate = 0.0;
        const Group2& groupParent = schedule.getGroup2(group.parent(), reportStepIdx);
        for (const std::string& groupName : groupParent.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);

            if ( (groupTmp.isProductionGroup() && !isInjector) ||
                 (groupTmp.isInjectionGroup() && isInjector) )
            {
                groupTotalGuideRate += guideRate->get(groupName, groupTarget);
            }
        }
        double groupGuideRate = guideRate->get(group.name(), groupTarget);
        return groupGuideRate / groupTotalGuideRate;
    }

    inline void accumulateGroupFractions(const std::string& groupName, const std::string& controlGroupName, const Schedule& schedule, const int reportStepIdx, const GuideRate* guideRate, const Group2::GuideRateTarget& groupTarget, const bool isInjector, double fraction) {

        const Group2& group = schedule.getGroup2(groupName, reportStepIdx);
        if (groupName != controlGroupName) {
            fraction *= groupFractionFromGuideRates(group, schedule, reportStepIdx, guideRate, groupTarget, isInjector);
            accumulateGroupFractions(group.parent(), controlGroupName, schedule, reportStepIdx, guideRate, groupTarget, isInjector, fraction);
        }

        return;
    }



    } // namespace wellGroupHelpers

}

#endif
