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
#include <opm/simulators/wells/WellAssemble.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellGroupControls.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceFluidSystem.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace Opm
{

template<class FluidSystem>
WellAssemble<FluidSystem>::
WellAssemble(const WellInterfaceFluidSystem<FluidSystem>& well)
    : well_(well)
{}

template<class FluidSystem>
template<class EvalWell>
void
WellAssemble<FluidSystem>::
assembleControlEqProd(const WellState& well_state,
                      const GroupState& group_state,
                      const Schedule& schedule,
                      const SummaryState& summaryState,
                      const Well::ProductionControls& controls,
                      const EvalWell& bhp,
                      const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                      const std::function<EvalWell()>& bhp_from_thp,
                      EvalWell& control_eq,
                      DeferredLogger& deferred_logger) const
{
    const auto current = well_state.well(well_.indexOfWell()).production_cmode;
    const auto& pu = well_.phaseUsage();
    const double efficiencyFactor = well_.wellEcl().getEfficiencyFactor();

    switch (current) {
    case Well::ProducerCMode::ORAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        const EvalWell rate = -rates[BlackoilPhases::Liquid];
        control_eq = rate - controls.oil_rate;
        break;
    }
    case Well::ProducerCMode::WRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        const EvalWell rate = -rates[BlackoilPhases::Aqua];
        control_eq = rate - controls.water_rate;
        break;
    }
    case Well::ProducerCMode::GRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
        const EvalWell rate = -rates[BlackoilPhases::Vapour];
        control_eq = rate - controls.gas_rate;
        break;
    }
    case Well::ProducerCMode::LRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        EvalWell rate = -rates[BlackoilPhases::Aqua] - rates[BlackoilPhases::Liquid];
        control_eq = rate - controls.liquid_rate;
        break;
    }
    case Well::ProducerCMode::CRAT: {
        OPM_DEFLOG_THROW(std::runtime_error,
                         "CRAT control not supported, well " + well_.name(),
                         deferred_logger);
    }
    case Well::ProducerCMode::RESV: {
        auto total_rate = rates[0]; // To get the correct type only.
        total_rate = 0.0;
        std::vector<double> convert_coeff(well_.numPhases(), 1.0);
        well_.rateConverter().calcCoeff(/*fipreg*/ 0, well_.pvtRegionIdx(), well_state.well(well_.indexOfWell()).surface_rates, convert_coeff);
        for (int phase = 0; phase < 3; ++phase) {
            if (pu.phase_used[phase]) {
                const int pos = pu.phase_pos[phase];
                total_rate -= rates[phase] * convert_coeff[pos]; // Note different indices.
            }
        }
        if (controls.prediction_mode) {
            control_eq = total_rate - controls.resv_rate;
        } else {
            std::vector<double> hrates(well_.numPhases(), 0.);
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                hrates[pu.phase_pos[Water]] = controls.water_rate;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                hrates[pu.phase_pos[Oil]] = controls.oil_rate;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                hrates[pu.phase_pos[Gas]] = controls.gas_rate;
            }
            std::vector<double> hrates_resv(well_.numPhases(), 0.);
            well_.rateConverter().calcReservoirVoidageRates(/*fipreg*/ 0, well_.pvtRegionIdx(), hrates, hrates_resv);
            double target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
            control_eq = total_rate - target;
        }
        break;
    }
    case Well::ProducerCMode::BHP: {
        control_eq = bhp - controls.bhp_limit;
        break;
    }
    case Well::ProducerCMode::THP: {
        control_eq = bhp - bhp_from_thp();
        break;
    }
    case Well::ProducerCMode::GRUP: {
        assert(well_.wellEcl().isAvailableForGroupControl());
        const auto& group = schedule.getGroup(well_.wellEcl().groupName(), well_.currentStep());
        // Annoying thing: the rates passed to this function are
        // always of size 3 and in canonical (for PhaseUsage)
        // order. This is what is needed for VFP calculations if
        // they are required (THP controlled well). But for the
        // group production control things we must pass only the
        // active phases' rates.
        std::vector<EvalWell> active_rates(pu.num_phases);
        for (int canonical_phase = 0; canonical_phase < 3; ++canonical_phase) {
            if (pu.phase_used[canonical_phase]) {
                active_rates[pu.phase_pos[canonical_phase]] = rates[canonical_phase];
            }
        }
        auto rCoeff = [this, &group_state](const RegionId id, const int region, const std::optional<std::string>& prod_gname, std::vector<double>& coeff)
        {
            if (prod_gname)
                well_.rateConverter().calcCoeff(id, region, group_state.production_rates(*prod_gname), coeff);
            else
                well_.rateConverter().calcCoeff(id, region, coeff);

        };
        WellGroupControls(well_).getGroupProductionControl(group, well_state,
                                                             group_state,
                                                             schedule,
                                                             summaryState,
                                                             bhp, active_rates,
                                                             rCoeff,
                                                             efficiencyFactor,
                                                             control_eq,
                                                             deferred_logger);
        break;
    }
    case Well::ProducerCMode::CMODE_UNDEFINED: {
        OPM_DEFLOG_THROW(std::runtime_error,
                         "Well control must be specified for well " + well_.name(),
                         deferred_logger);
    }
    case Well::ProducerCMode::NONE: {
        OPM_DEFLOG_THROW(std::runtime_error,
                         "Well control must be specified for well " + well_.name(),
                         deferred_logger);
    }
    }
}

template<class FluidSystem>
template<class EvalWell>
void
WellAssemble<FluidSystem>::
assembleControlEqInj(const WellState& well_state,
                     const GroupState& group_state,
                     const Schedule& schedule,
                     const SummaryState& summaryState,
                     const Well::InjectionControls& controls,
                     const EvalWell& bhp,
                     const EvalWell& injection_rate,
                     const std::function<EvalWell()>& bhp_from_thp,
                     EvalWell& control_eq,
                     DeferredLogger& deferred_logger) const
{
    auto current = well_state.well(well_.indexOfWell()).injection_cmode;
    const InjectorType injectorType = controls.injector_type;
    const auto& pu = well_.phaseUsage();
    const double efficiencyFactor = well_.wellEcl().getEfficiencyFactor();

    switch (current) {
    case Well::InjectorCMode::RATE: {
        control_eq = injection_rate - controls.surface_rate;
        break;
    }
    case Well::InjectorCMode::RESV: {
        std::vector<double> convert_coeff(well_.numPhases(), 1.0);
        well_.rateConverter().calcInjCoeff(/*fipreg*/ 0, well_.pvtRegionIdx(), convert_coeff);

        double coeff;

        switch (injectorType) {
        case InjectorType::WATER: {
            coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Aqua]];
            break;
        }
        case InjectorType::OIL: {
            coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Liquid]];
            break;
        }
        case InjectorType::GAS: {
            coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Vapour]];
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injectors " + well_.wellEcl().name());
        }

        control_eq = coeff * injection_rate - controls.reservoir_rate;
        break;
    }
    case Well::InjectorCMode::THP: {
        control_eq = bhp - bhp_from_thp();
        break;
    }
    case Well::InjectorCMode::BHP: {
        control_eq = bhp - controls.bhp_limit;
        break;
    }
    case Well::InjectorCMode::GRUP: {
        assert(well_.wellEcl().isAvailableForGroupControl());
        const auto& group = schedule.getGroup(well_.wellEcl().groupName(), well_.currentStep());
        auto rCoeff = [this, &group_state](const RegionId id, const int region, const std::optional<std::string>& prod_gname, std::vector<double>& coeff)
        {
            if(prod_gname) {
                well_.rateConverter().calcCoeff(id, region, group_state.production_rates(*prod_gname), coeff);
            } else {
                well_.rateConverter().calcInjCoeff(id, region, coeff);
            }
        };
        WellGroupControls(well_).getGroupInjectionControl(group,
                                                            well_state,
                                                            group_state,
                                                            schedule,
                                                            summaryState,
                                                            injectorType,
                                                            bhp,
                                                            injection_rate,
                                                            rCoeff,
                                                            efficiencyFactor,
                                                            control_eq,
                                                            deferred_logger);
        break;
    }
    case Well::InjectorCMode::CMODE_UNDEFINED: {
        OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + well_.name(), deferred_logger);
    }
    }
}

#define INSTANCE_METHODS(A,...) \
template void WellAssemble<A>:: \
assembleControlEqProd<__VA_ARGS__>(const WellState&, \
                                   const GroupState&, \
                                   const Schedule&, \
                                   const SummaryState&, \
                                   const Well::ProductionControls&, \
                                   const __VA_ARGS__&, \
                                   const std::vector<__VA_ARGS__>&, \
                                   const std::function<__VA_ARGS__()>&, \
                                   __VA_ARGS__&, \
                                   DeferredLogger&) const; \
template void WellAssemble<A>:: \
assembleControlEqInj<__VA_ARGS__>(const WellState&, \
                                  const GroupState&, \
                                  const Schedule&, \
                                  const SummaryState&, \
                                  const Well::InjectionControls&, \
                                  const __VA_ARGS__&, \
                                  const __VA_ARGS__&, \
                                  const std::function<__VA_ARGS__()>&, \
                                  __VA_ARGS__&, \
                                  DeferredLogger&) const;

using FluidSys = BlackOilFluidSystem<double, BlackOilDefaultIndexTraits>;

template class WellAssemble<FluidSys>;

INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,3,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,4,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,5,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,6,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,7,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,8,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,9,0u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,4u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,5u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,6u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,7u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,8u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,9u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,10u>)
INSTANCE_METHODS(FluidSys, DenseAd::Evaluation<double,-1,11u>)
} // namespace Opm
