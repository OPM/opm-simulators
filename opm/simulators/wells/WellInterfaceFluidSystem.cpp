/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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

#include <config.h>
#include <opm/simulators/wells/WellInterfaceFluidSystem.hpp>

#include <opm/grid/utility/RegionMapping.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
#include <opm/simulators/wells/WellConstraints.hpp>
#include <opm/simulators/wells/WellGroupConstraints.hpp>
#include <opm/simulators/wells/WellGroupControls.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>
#include <cmath>

namespace Opm
{

template<class FluidSystem>
WellInterfaceFluidSystem<FluidSystem>::
WellInterfaceFluidSystem(const Well& well,
                         const ParallelWellInfo& parallel_well_info,
                         const int time_step,
                         const RateConverterType& rate_converter,
                         const int pvtRegionIdx,
                         const int num_components,
                         const int num_phases,
                         const int index_of_well,
                         const std::vector<PerforationData>& perf_data)
    : WellInterfaceGeneric(well, parallel_well_info, time_step,
                           pvtRegionIdx, num_components, num_phases,
                           index_of_well, perf_data)
    , rateConverter_(rate_converter)
{
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
calculateReservoirRates(SingleWellState& ws) const
{
    const int fipreg = 0; // not considering the region for now
    const int np = number_of_phases_;

    std::vector<double> surface_rates(np, 0.0);
    for (int p = 0; p < np; ++p) {
        surface_rates[p] = ws.surface_rates[p];
    }

    std::vector<double> voidage_rates(np, 0.0);
    rateConverter_.calcReservoirVoidageRates(fipreg, pvtRegionIdx_, surface_rates, voidage_rates);
    ws.reservoir_rates = voidage_rates;
}

template <typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkIndividualConstraints(SingleWellState& ws,
                           const SummaryState& summaryState,
                           DeferredLogger& deferred_logger) const
{
    auto rRates = [this](const int fipreg,
                         const int pvtRegion,
                         const std::vector<double>& surface_rates,
                         std::vector<double>& voidage_rates)
    {
        return rateConverter_.calcReservoirVoidageRates(fipreg, pvtRegion,
                                                        surface_rates, voidage_rates);
    };

    return WellConstraints(*this).
            checkIndividualConstraints(ws, summaryState, rRates,
                                       this->operability_status_.thp_limit_violated_but_not_switched,
                                       deferred_logger);
}

template <typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkGroupConstraints(WellState& well_state,
                      const GroupState& group_state,
                      const Schedule& schedule,
                      const SummaryState& summaryState,
                      DeferredLogger& deferred_logger) const
{
    auto rCoeff = [this](const int id, const int region, std::vector<double>& coeff)
    {
        this->rateConverter().calcCoeff(id, region, coeff);
    };

    return WellGroupConstraints(*this).checkGroupConstraints(well_state, group_state,
                                                             schedule, summaryState,
                                                             rCoeff, deferred_logger);
}

template <typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkConstraints(WellState& well_state,
                 const GroupState& group_state,
                 const Schedule& schedule,
                 const SummaryState& summaryState,
                 DeferredLogger& deferred_logger) const
{
    const bool ind_broken = checkIndividualConstraints(well_state.well(this->index_of_well_), summaryState, deferred_logger);
    if (ind_broken) {
        return true;
    } else {
        return checkGroupConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
    }
}

template<typename FluidSystem>
int
WellInterfaceFluidSystem<FluidSystem>::
flowPhaseToEbosPhaseIdx(const int phaseIdx) const
{
    const auto& pu = this->phaseUsage();
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
        return FluidSystem::waterPhaseIdx;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
        return FluidSystem::oilPhaseIdx;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
        return FluidSystem::gasPhaseIdx;

    // for other phases return the index
    return phaseIdx;
}

template<typename FluidSystem>
std::optional<double>
WellInterfaceFluidSystem<FluidSystem>::
getGroupInjectionTargetRate(const Group& group,
                            const WellState& well_state,
                            const GroupState& group_state,
                            const Schedule& schedule,
                            const SummaryState& summaryState,
                            const InjectorType& injectorType,
                            double efficiencyFactor,
                            DeferredLogger& deferred_logger) const
{
    auto rCoeff = [this](const int id, const int region, std::vector<double>& coeff)
    {
        this->rateConverter().calcCoeff(id, region, coeff);
    };

    return WellGroupControls(*this).getGroupInjectionTargetRate(group, well_state,
                                                                group_state, schedule,
                                                                summaryState, injectorType,
                                                                rCoeff, efficiencyFactor,
                                                                deferred_logger);
}

template<typename FluidSystem>
double
WellInterfaceFluidSystem<FluidSystem>::
getGroupProductionTargetRate(const Group& group,
                          const WellState& well_state,
                          const GroupState& group_state,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          double efficiencyFactor) const
{
    auto rCoeff = [this](const int id, const int region, std::vector<double>& coeff)
    {
        this->rateConverter().calcCoeff(id, region, coeff);
    };

    return WellGroupControls(*this).getGroupProductionTargetRate(group, well_state,
                                                                 group_state, schedule,
                                                                 summaryState,
                                                                 rCoeff, efficiencyFactor);
}

template class WellInterfaceFluidSystem<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>>;

} // namespace Opm
