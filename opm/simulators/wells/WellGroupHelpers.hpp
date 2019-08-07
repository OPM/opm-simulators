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


    inline double sumWellRates(const Group2& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos, const bool injector) {

        double rate = 0.0;
        for (const std::string& groupName : group.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);
            rate += groupTmp.getGroupEfficiencyFactor()*sumWellRates(groupTmp, schedule, wellState, reportStepIdx, phasePos, injector);
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

            double factor = wellEcl.getEfficiencyFactor();
            const auto wellrate_index = well_index * wellState.numPhases();
            if (injector)
                rate += factor * wellState.wellRates()[ wellrate_index + phasePos];
            else
                rate -= factor * wellState.wellRates()[ wellrate_index + phasePos];
        }
        return rate;
    }

    inline double sumWellResRates(const Group2& group, const Schedule& schedule, const WellStateFullyImplicitBlackoil& wellState, const int reportStepIdx, const int phasePos, const bool injector) {

        double rate = 0.0;
        for (const std::string& groupName : group.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);
            rate += groupTmp.getGroupEfficiencyFactor() * sumWellResRates(groupTmp, schedule, wellState, reportStepIdx, phasePos, injector);
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

            double factor = wellEcl.getEfficiencyFactor();

            const auto wellrate_index = well_index * wellState.numPhases();
            if (injector)
                rate += factor * wellState.wellReservoirRates()[ wellrate_index + phasePos];
            else
                rate -= factor * wellState.wellReservoirRates()[ wellrate_index + phasePos];
        }
        return rate;
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
            if (!wellEcl.isAvailableForGroupControl())
                continue;
            if (wellEcl.isProducer() && !injector)
                wellState.currentProductionControls()[well_index] = Well2::ProducerCMode::GRUP;

            if (wellEcl.isInjector() && injector)
                wellState.currentInjectionControls()[well_index] = Well2::InjectorCMode::GRUP;
        }
    }

    inline double computeGuideRate(const WellStateFullyImplicitBlackoil& wellState, const Well2& well, const int phasePos )
    {                
        const auto& end = wellState.wellMap().end();
        const auto& it = wellState.wellMap().find( well.name() );
        if (it == end)  // the well is not found
            return 0.0;

        int well_index = it->second[0];

         Opm::Well2::GuideRateTarget GuideRatePhase = well.getGuideRatePhase();
#warning handle the case when phase != GuideRatePhase
        if (GuideRatePhase == Opm::Well2::GuideRateTarget::UNDEFINED) {
            const auto wellrate_index = well_index * wellState.numPhases();
            double guideRate = wellState.wellPotentials()[wellrate_index + phasePos];
            //if (guideRate == 0)
            //    guideRate = 1; // if the well potential is not defined set it to 1
#warning Add support for GUIDERAT
            return guideRate*well.getGuideRateScalingFactor();
        }
        return well.getGuideRate()*well.getGuideRateScalingFactor();
    }

    inline void computeTotalGuideRate(const Group2& group, const WellStateFullyImplicitBlackoil& wellState, const Schedule& schedule, const int reportStepIdx, const int phasePos, const bool isInjector, double& groupTotalGuideRate, double& groupTargetReduction )
    {

        for (const std::string& groupName : group.groups()) {
            const Group2& groupTmp = schedule.getGroup2(groupName, reportStepIdx);
            computeTotalGuideRate(groupTmp, wellState, schedule, reportStepIdx, phasePos, isInjector, groupTotalGuideRate, groupTargetReduction);
        }
#warning return groupGuideRate if set, else sum it from below
        for (const std::string& wellName : group.wells()) {
            const auto& wellTmp = schedule.getWell2(wellName, reportStepIdx);

            if (wellTmp.isProducer() && isInjector)
                 continue;

            if (wellTmp.isInjector() && !isInjector)
                 continue;

            const auto& end = wellState.wellMap().end();
                const auto& it = wellState.wellMap().find( wellName );
                if (it == end)  // the well is not found
                    continue;

            int well_index = it->second[0];
            const auto wellrate_index = well_index * wellState.numPhases();
            // only include wells under group control
            if (isInjector) {
                if (wellState.currentInjectionControls()[well_index] == Well2::InjectorCMode::GRUP)
                    groupTotalGuideRate += computeGuideRate(wellState, wellTmp, phasePos);
                else
                    groupTargetReduction += wellState.wellRates()[wellrate_index + phasePos];
            } else {
                if (wellState.currentProductionControls()[well_index] ==  Well2::ProducerCMode::GRUP)
                    groupTotalGuideRate += computeGuideRate(wellState, wellTmp, phasePos);
                else
                    groupTargetReduction -= wellState.wellRates()[wellrate_index + phasePos];
            }
        }
    }


    } // namespace wellGroupHelpers

}

#endif
