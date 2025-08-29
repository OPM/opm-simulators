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
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>

#include <cassert>
#include <stdexcept>

namespace Opm
{

template<typename FluidSystem>
WellAssemble<FluidSystem>::
WellAssemble(const WellInterfaceFluidSystem<FluidSystem>& well)
    : well_(well)
{}

template<typename FluidSystem>
template<class EvalWell>
void
WellAssemble<FluidSystem>::
assembleControlEqProd(const WellState<Scalar, IndexTraits>& well_state,
                      const GroupState<Scalar>& group_state,
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
    const Scalar efficiencyFactor = well_.wellEcl().getEfficiencyFactor() *
                                    well_state[well_.name()].efficiency_scaling_factor;

    switch (current) {
    case Well::ProducerCMode::ORAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        const EvalWell rate = -rates[FluidSystem::oilPhaseIdx];
        control_eq = rate - controls.oil_rate;
        break;
    }
    case Well::ProducerCMode::WRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        const EvalWell rate = -rates[FluidSystem::waterPhaseIdx];
        control_eq = rate - controls.water_rate;
        break;
    }
    case Well::ProducerCMode::GRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
        const EvalWell rate = -rates[FluidSystem::gasPhaseIdx];
        control_eq = rate - controls.gas_rate;
        break;
    }
    case Well::ProducerCMode::LRAT: {
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        EvalWell rate = -rates[FluidSystem::waterPhaseIdx] - rates[FluidSystem::oilPhaseIdx];
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
        std::vector<Scalar> convert_coeff(well_.numPhases(), 1.0);
        well_.rateConverter().calcCoeff(/*fipreg*/ 0, well_.pvtRegionIdx(), well_state.well(well_.indexOfWell()).surface_rates, convert_coeff);
        for (int phase = 0; phase < 3; ++phase) {
            if (FluidSystem::phaseIsActive(phase)) {
                const int pos = FluidSystem::canonicalToActivePhaseIdx(phase);
                total_rate -= rates[phase] * convert_coeff[pos]; // Note different indices.
            }
        }
        if (controls.prediction_mode) {
            control_eq = total_rate - controls.resv_rate;
        } else {
            std::vector<Scalar> hrates(well_.numPhases(), 0.);
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const int water_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
                hrates[water_pos] = controls.water_rate;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const int oil_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
                hrates[oil_pos] = controls.oil_rate;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const int gas_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
                hrates[gas_pos] = controls.gas_rate;
            }
            std::vector<Scalar> hrates_resv(well_.numPhases(), 0.);
            well_.rateConverter().calcReservoirVoidageRates(/*fipreg*/ 0, well_.pvtRegionIdx(), hrates, hrates_resv);
            Scalar target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
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
        const auto& pu = well_.phaseUsage();
        std::vector<EvalWell> active_rates(pu.numActivePhases());
        for (int canonical_phase = 0; canonical_phase < 3; ++canonical_phase) {
            if (FluidSystem::phaseIsActive(canonical_phase)) {
                const int phase_pos = FluidSystem::canonicalToActivePhaseIdx(canonical_phase);
                active_rates[phase_pos] = rates[canonical_phase];
            }
        }
        auto rCoeff = [this, &group_state](const RegionId id, const int region,
                                           const std::optional<std::string>& prod_gname,
                                           std::vector<Scalar>& coeff)
        {
            if (prod_gname)
                well_.rateConverter().calcCoeff(id, region, group_state.production_rates(*prod_gname), coeff);
            else
                well_.rateConverter().calcCoeff(id, region, coeff);

        };

        WellGroupControls(well_).getGroupProductionControl(group, 
                                                           well_state,
                                                           group_state,
                                                           schedule,
                                                           summaryState,
                                                           bhp, 
                                                           active_rates,
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

template<typename FluidSystem>
template<class EvalWell>
void
WellAssemble<FluidSystem>::
assembleControlEqInj(const WellState<Scalar, IndexTraits>& well_state,
                     const GroupState<Scalar>& group_state,
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
    const Scalar efficiencyFactor = well_.wellEcl().getEfficiencyFactor() *
                                    well_state[well_.name()].efficiency_scaling_factor;

    switch (current) {
    case Well::InjectorCMode::RATE: {
        control_eq = injection_rate - controls.surface_rate;
        break;
    }
    case Well::InjectorCMode::RESV: {
        std::vector<Scalar> convert_coeff(well_.numPhases(), 1.0);
        well_.rateConverter().calcInjCoeff(/*fipreg*/ 0, well_.pvtRegionIdx(), convert_coeff);

        Scalar coeff;

        switch (injectorType) {
        case InjectorType::WATER: {
            const int phase_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
            coeff = convert_coeff[phase_pos];
            break;
        }
        case InjectorType::OIL: {
            const int phase_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
            coeff = convert_coeff[phase_pos];
            break;
        }
        case InjectorType::GAS: {
            const int phase_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
            coeff = convert_coeff[phase_pos];
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
        auto rCoeff = [this, &group_state](const RegionId id, const int region,
                                           const std::optional<std::string>& prod_gname,
                                           std::vector<Scalar>& coeff)
        {
            if(prod_gname) {
                well_.rateConverter().calcCoeff(id, region,
                                                group_state.production_rates(*prod_gname), coeff);
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

#define INSTANTIATE_METHODS(A,...)                                        \
template void WellAssemble<A>::                                           \
assembleControlEqProd<__VA_ARGS__>(const WellState<typename A::Scalar, BlackOilDefaultFluidSystemIndices>&,  \
                                   const GroupState<typename A::Scalar>&, \
                                   const Schedule&,                       \
                                   const SummaryState&,                   \
                                   const Well::ProductionControls&,       \
                                   const __VA_ARGS__&,                    \
                                   const std::vector<__VA_ARGS__>&,       \
                                   const std::function<__VA_ARGS__()>&,   \
                                   __VA_ARGS__&,                          \
                                   DeferredLogger&) const;                \
template void WellAssemble<A>::                                           \
assembleControlEqInj<__VA_ARGS__>(const WellState<typename A::Scalar, BlackOilDefaultFluidSystemIndices>&,   \
                                  const GroupState<typename A::Scalar>&,  \
                                  const Schedule&,                        \
                                  const SummaryState&,                    \
                                  const Well::InjectionControls&,         \
                                  const __VA_ARGS__&,                     \
                                  const __VA_ARGS__&,                     \
                                  const std::function<__VA_ARGS__()>&,    \
                                  __VA_ARGS__&,                           \
                                  DeferredLogger&) const;

template<class Scalar>
using FS = BlackOilFluidSystem<Scalar, BlackOilDefaultFluidSystemIndices>;

#define INSTANTIATE_TYPE(T)                                   \
    template class WellAssemble<FS<T>>;                       \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,3,0u>)   \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,4,0u>)   \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,5,0u>)   \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,6,0u>)   \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,7,0u>)   \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,8,0u>)   \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,9,0u>)   \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,4u>)  \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,5u>)  \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,6u>)  \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,7u>)  \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,8u>)  \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,9u>)  \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,10u>) \
    INSTANTIATE_METHODS(FS<T>, DenseAd::Evaluation<T,-1,11u>)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
