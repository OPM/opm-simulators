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

template <typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
calculateReservoirRates(SingleWellState& ws) const
{
    const int np = number_of_phases_;
    const auto& pu = this->phaseUsage();
    // Calculate reservoir rates from average pressure and temperature
    if ( !pu.has_energy || this->wellEcl().isProducer()) {
        const int fipreg = 0; // not considering the region for now
        this->rateConverter_
            .calcReservoirVoidageRates(fipreg,
                                    this->pvtRegionIdx_,
                                    ws.surface_rates,
                                    ws.reservoir_rates);

        // Compute total connection reservoir rate CVPR/CVIR
        auto& perf_data = ws.perf_data;
        const auto num_perf_well = perf_data.size();
        const auto& surf_perf_rates = perf_data.phase_rates;
        for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
            const auto surface_rates_perf = std::vector<double>
                { surf_perf_rates.begin() + (i + 0)*np ,
                surf_perf_rates.begin() + (i + 1)*np };

            std::vector<double> voidage_rates_perf(np, 0.0);
            this->rateConverter_
                .calcReservoirVoidageRates(fipreg,
                                        this->pvtRegionIdx_,
                                        surface_rates_perf,
                                        voidage_rates_perf);

            perf_data.rates[i] =
                std::accumulate(voidage_rates_perf.begin(),
                                voidage_rates_perf.end(), 0.0);
        }
        return;
    }
    // For injectors in a thermal case we convert using the well bhp and temperature
    // Assume pure phases in the injector
    const auto saltConc = 0.0;
    auto rsMax = 0.0;
    auto rvMax = 0.0;
    auto rswMax = 0.0;
    auto rvwMax = 0.0;
    this->rateConverter_
        .calcReservoirVoidageRates(this->pvtRegionIdx_,
                                    ws.bhp,
                                    rsMax,
                                    rvMax,
                                    rswMax,
                                    rvwMax,
                                    ws.temperature,
                                    saltConc,
                                    ws.surface_rates,
                                    ws.reservoir_rates);


    // Compute total connection reservoir rate CVIR
    auto& perf_data = ws.perf_data;
    const auto num_perf_well = perf_data.size();
    const auto& surf_perf_rates = perf_data.phase_rates;
    for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
        const auto surface_rates_perf = std::vector<double>
            { surf_perf_rates.begin() + (i + 0)*np ,
              surf_perf_rates.begin() + (i + 1)*np };

        const auto pressure = perf_data.pressure[i];
        // Calculate other per-phase dynamic quantities.
        const auto temperature = ws.temperature; // Assume same  temperature in the well
        std::vector<double> voidage_rates_perf(np, 0.0);
        this->rateConverter_
            .calcReservoirVoidageRates(this->pvtRegionIdx_,
                                       pressure,
                                       rsMax,
                                       rvMax,
                                       rswMax, // Rsw
                                       rvwMax, // Rvw
                                       temperature,
                                       saltConc,
                                       surface_rates_perf,
                                       voidage_rates_perf);

        perf_data.rates[i] =
            std::accumulate(voidage_rates_perf.begin(),
                            voidage_rates_perf.end(), 0.0);
    }
}

template <typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkIndividualConstraints(SingleWellState& ws,
                           const SummaryState& summaryState,
                           DeferredLogger& deferred_logger,
                           const std::optional<Well::InjectionControls>& inj_controls,
                           const std::optional<Well::ProductionControls>& prod_controls) const
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
                                       deferred_logger, inj_controls, prod_controls);
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

    if (!this->wellEcl().isAvailableForGroupControl())
        return false;

    auto rCoeff = [this, &group_state](const RegionId id, const int region, const std::optional<std::string>& prod_gname, std::vector<double>& coeff)
    {
        if (prod_gname)
            this->rateConverter().calcCoeff(id, region, group_state.production_rates(*prod_gname), coeff);
        else
            this->rateConverter().calcInjCoeff(id, region, coeff);

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
    auto rCoeff = [this, &group_state](const RegionId id, const int region, const std::optional<std::string>& prod_gname, std::vector<double>& coeff)
    {
        if (prod_gname)
            this->rateConverter().calcCoeff(id, region, group_state.production_rates(*prod_gname), coeff);
        else
            this->rateConverter().calcInjCoeff(id, region, coeff);

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
                          double efficiencyFactor,
                          DeferredLogger& deferred_logger) const
{
    auto rCoeff = [this, &group_state](const RegionId id, const int region, const std::optional<std::string>& prod_gname, std::vector<double>& coeff)
    {
        if (prod_gname)
            this->rateConverter().calcCoeff(id, region, group_state.production_rates(*prod_gname), coeff);
        else
            this->rateConverter().calcInjCoeff(id, region, coeff);

    };

    return WellGroupControls(*this).getGroupProductionTargetRate(group, well_state,
                                                                 group_state, schedule,
                                                                 summaryState,
                                                                 rCoeff, efficiencyFactor,
                                                                 deferred_logger);
}

template class WellInterfaceFluidSystem<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>>;

} // namespace Opm
