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

namespace Opm
{

template<typename FluidSystem>
WellInterfaceFluidSystem<FluidSystem>::
WellInterfaceFluidSystem(const Well& well,
                         const ParallelWellInfo<Scalar>& parallel_well_info,
                         const int time_step,
                         const ModelParameters& param,
                         const RateConverterType& rate_converter,
                         const int pvtRegionIdx,
                         const int num_conservation_quantities,
                         const int num_phases,
                         const int index_of_well,
                         const std::vector<PerforationData<Scalar>>& perf_data)
    : WellInterfaceGeneric<Scalar, IndexTraits>(well, parallel_well_info, time_step, param,
                                   pvtRegionIdx, num_conservation_quantities, num_phases,
                                   index_of_well, FluidSystem::phaseUsage(), perf_data)
    , rateConverter_(rate_converter)
{
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
calculateReservoirRates(const bool co2store, SingleWellState<Scalar, IndexTraits>& ws) const
{
    const int np = this->number_of_phases_;
    const auto& pu = this->phaseUsage();
    // Calculate reservoir rates from average pressure and temperature
    if ( !(co2store || pu.hasEnergy()) || this->wellEcl().isProducer()) {
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
            const auto surface_rates_perf = std::vector<Scalar>
                { surf_perf_rates.begin() + (i + 0)*np ,
                surf_perf_rates.begin() + (i + 1)*np };

            std::vector<Scalar> voidage_rates_perf(np, 0.0);
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
    // For injectors in a co2 storage case or a thermal case
    // we convert using the well bhp and temperature
    // Assume pure phases in the injector
    const Scalar saltConc = 0.0;
    Scalar rsMax = 0.0;
    Scalar rvMax = 0.0;
    Scalar rswMax = 0.0;
    Scalar rvwMax = 0.0;
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
        const auto surface_rates_perf = std::vector<Scalar>
            { surf_perf_rates.begin() + (i + 0)*np ,
              surf_perf_rates.begin() + (i + 1)*np };

        const auto pressure = perf_data.pressure[i];
        // Calculate other per-phase dynamic quantities.
        const auto temperature = ws.temperature; // Assume same  temperature in the well
        std::vector<Scalar> voidage_rates_perf(np, 0.0);
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

template<typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkIndividualConstraints(SingleWellState<Scalar, IndexTraits>& ws,
                           const SummaryState& summaryState,
                           DeferredLogger& deferred_logger,
                           const std::optional<Well::InjectionControls>& inj_controls,
                           const std::optional<Well::ProductionControls>& prod_controls) const
{
    auto rRates = [this](const int fipreg,
                         const int pvtRegion,
                         const std::vector<Scalar>& surface_rates,
                         std::vector<Scalar>& voidage_rates)
    {
        return rateConverter_.calcReservoirVoidageRates(fipreg, pvtRegion,
                                                        surface_rates, voidage_rates);
    };

    return WellConstraints(*this).
            checkIndividualConstraints(ws, summaryState, rRates,
                                       this->operability_status_.thp_limit_violated_but_not_switched,
                                       deferred_logger, inj_controls, prod_controls);
}

template<typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkGroupConstraints(WellState<Scalar, IndexTraits>& well_state,
                      const GroupState<Scalar>& group_state,
                      const Schedule& schedule,
                      const SummaryState& summaryState,
                      const bool check_guide_rate,
                      DeferredLogger& deferred_logger) const
{
    if (!this->wellEcl().isAvailableForGroupControl())
        return false;

    auto rCoeff = [this, &group_state](const RegionId id,
                                       const int region,
                                       const std::optional<std::string>& prod_gname,
                                       std::vector<Scalar>& coeff)
    {
        if (prod_gname)
            this->rateConverter().calcCoeff(id, region,
                                            group_state.production_rates(*prod_gname), coeff);
        else
            this->rateConverter().calcInjCoeff(id, region, coeff);
    };

    return WellGroupConstraints(*this).checkGroupConstraints(well_state, group_state,
                                                             schedule, summaryState,
                                                             rCoeff, check_guide_rate, deferred_logger);
}

template<typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkConstraints(WellState<Scalar, IndexTraits>& well_state,
                 const GroupState<Scalar>& group_state,
                 const Schedule& schedule,
                 const SummaryState& summaryState,
                 DeferredLogger& deferred_logger) const
{
    const bool ind_broken = checkIndividualConstraints(well_state.well(this->index_of_well_),
                                                       summaryState, deferred_logger);
    if (ind_broken) {
        return true;
    } else {
        return checkGroupConstraints(well_state, group_state, schedule,
                                     summaryState, true, deferred_logger);
    }
}

template<typename FluidSystem>
std::optional<typename FluidSystem::Scalar>
WellInterfaceFluidSystem<FluidSystem>::
getGroupInjectionTargetRate(const Group& group,
                            const WellState<Scalar, IndexTraits>& well_state,
                            const GroupState<Scalar>& group_state,
                            const Schedule& schedule,
                            const SummaryState& summaryState,
                            const InjectorType& injectorType,
                            Scalar efficiencyFactor,
                            DeferredLogger& deferred_logger) const
{
    auto rCoeff = [this, &group_state](const RegionId id, const int region,
                                       const std::optional<std::string>& prod_gname,
                                       std::vector<Scalar>& coeff)
    {
        if (prod_gname)
            this->rateConverter().calcCoeff(id, region,
                                            group_state.production_rates(*prod_gname), coeff);
        else
            this->rateConverter().calcInjCoeff(id, region, coeff);
    };

    return WellGroupControls(*this).getGroupInjectionTargetRate(group,
                                                                well_state,
                                                                group_state,
                                                                schedule,
                                                                summaryState,
                                                                injectorType,
                                                                rCoeff,
                                                                efficiencyFactor,
                                                                deferred_logger);
}

template<typename FluidSystem>
typename FluidSystem::Scalar
WellInterfaceFluidSystem<FluidSystem>::
getGroupProductionTargetRate(const Group& group,
                             const WellState<Scalar, IndexTraits>& well_state,
                             const GroupState<Scalar>& group_state,
                             const Schedule& schedule,
                             const SummaryState& summaryState,
                             Scalar efficiencyFactor,
                             DeferredLogger& deferred_logger) const
{
    auto rCoeff = [this, &group_state](const RegionId id, const int region,
                                       const std::optional<std::string>& prod_gname,
                                       std::vector<Scalar>& coeff)
    {
        if (prod_gname)
            this->rateConverter().calcCoeff(id, region,
                                            group_state.production_rates(*prod_gname), coeff);
        else
            this->rateConverter().calcInjCoeff(id, region, coeff);
    };

    return WellGroupControls(*this).getGroupProductionTargetRate(group,
                                                                 well_state,
                                                                 group_state,
                                                                 schedule,
                                                                 summaryState,
                                                                 rCoeff,
                                                                 efficiencyFactor,
                                                                 deferred_logger);
}

template<typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
zeroGroupRateTarget(const SummaryState& summary_state,
                    const Schedule& schedule,
                    const WellState<Scalar, IndexTraits>& well_state,
                    const GroupState<Scalar>& group_state,
                    DeferredLogger& deferred_logger) const
{
    const auto& well = this->well_ecl_;
    const auto& group = schedule.getGroup(well.groupName(), this->currentStep());
    const Scalar efficiencyFactor = well.getEfficiencyFactor() *
                                    well_state[well.name()].efficiency_scaling_factor;
    if (this->isInjector()) {
        // Check injector under group control
        const auto& controls = well.injectionControls(summary_state);
        const std::optional<Scalar> target =
            this->getGroupInjectionTargetRate(group, well_state,
                                              group_state, schedule,
                                              summary_state, controls.injector_type,
                                              efficiencyFactor, deferred_logger);
        if (target.has_value()) {
            return target.value() == 0.0;
        } else {
            return false;
        }
    } else {
        // Check producer under group control
        const Scalar scale =
            this->getGroupProductionTargetRate(group, well_state,
                                               group_state, schedule,
                                               summary_state, efficiencyFactor,
                                               deferred_logger);
        return scale == 0.0;
    }
}

template<class Scalar>
using FS = BlackOilFluidSystem<Scalar, BlackOilDefaultFluidSystemIndices>;

template class WellInterfaceFluidSystem<FS<double>>;

#if FLOW_INSTANTIATE_FLOAT
template class WellInterfaceFluidSystem<FS<float>>;
#endif

} // namespace Opm
