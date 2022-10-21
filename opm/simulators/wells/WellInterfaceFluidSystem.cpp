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

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>
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
Well::ProducerCMode
WellInterfaceFluidSystem<FluidSystem>::
activeProductionConstraint(const SingleWellState& ws,
                           const SummaryState& summaryState,
                           DeferredLogger& deferred_logger) const
{
    const PhaseUsage& pu = this->phaseUsage();
    const auto controls = this->well_ecl_.productionControls(summaryState);
    const auto currentControl = ws.production_cmode;

    if (controls.hasControl(Well::ProducerCMode::BHP) && currentControl != Well::ProducerCMode::BHP) {
        const double bhp_limit = controls.bhp_limit;
        double current_bhp = ws.bhp;
        if (bhp_limit > current_bhp)
            return Well::ProducerCMode::BHP;
    }

    if (controls.hasControl(Well::ProducerCMode::ORAT) && currentControl != Well::ProducerCMode::ORAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]];
        if (controls.oil_rate < current_rate)
            return Well::ProducerCMode::ORAT;
    }

    if (controls.hasControl(Well::ProducerCMode::WRAT) && currentControl != Well::ProducerCMode::WRAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];
        if (controls.water_rate < current_rate)
            return Well::ProducerCMode::WRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::GRAT) && currentControl != Well::ProducerCMode::GRAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Vapour]];
        if (controls.gas_rate < current_rate)
            return Well::ProducerCMode::GRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::LRAT) && currentControl != Well::ProducerCMode::LRAT) {
        double current_rate = -ws.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]];
        current_rate -= ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];

        bool skip = false;
        if (controls.liquid_rate == controls.oil_rate) {
            const double current_water_rate = ws.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]];
            if (std::abs(current_water_rate) < 1e-12) {
                skip = true;
                deferred_logger.debug("LRAT_ORAT_WELL", "Well " + this->name() + " The LRAT target is equal the ORAT target and the water rate is zero, skip checking LRAT");
            }
        }
        if (!skip && controls.liquid_rate < current_rate)
            return Well::ProducerCMode::LRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::RESV) && currentControl != Well::ProducerCMode::RESV) {
        double current_rate = 0.0;
        if (pu.phase_used[BlackoilPhases::Aqua])
            current_rate -= ws.reservoir_rates[pu.phase_pos[BlackoilPhases::Aqua]];

        if (pu.phase_used[BlackoilPhases::Liquid])
            current_rate -= ws.reservoir_rates[pu.phase_pos[BlackoilPhases::Liquid]];

        if (pu.phase_used[BlackoilPhases::Vapour])
            current_rate -= ws.reservoir_rates[pu.phase_pos[BlackoilPhases::Vapour]];

        if (controls.prediction_mode && controls.resv_rate < current_rate)
            return Well::ProducerCMode::RESV;

        if (!controls.prediction_mode) {
            const int fipreg = 0; // not considering the region for now
            const int np = number_of_phases_;

            std::vector<double> surface_rates(np, 0.0);
            if (pu.phase_used[BlackoilPhases::Aqua])
                surface_rates[pu.phase_pos[BlackoilPhases::Aqua]] = controls.water_rate;
            if (pu.phase_used[BlackoilPhases::Liquid])
                surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = controls.oil_rate;
            if (pu.phase_used[BlackoilPhases::Vapour])
                surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = controls.gas_rate;

            std::vector<double> voidage_rates(np, 0.0);
            rateConverter_.calcReservoirVoidageRates(fipreg, pvtRegionIdx_, surface_rates, voidage_rates);

            double resv_rate = 0.0;
            for (int p = 0; p < np; ++p)
                resv_rate += voidage_rates[p];

            if (resv_rate < current_rate)
                return Well::ProducerCMode::RESV;
        }
    }

    if (controls.hasControl(Well::ProducerCMode::THP) && currentControl != Well::ProducerCMode::THP) {
        const auto& thp = getTHPConstraint(summaryState);
        double current_thp = ws.thp;
        if (thp > current_thp && !ws.trivial_target) {
            // If WVFPEXP item 4 is set to YES1 or YES2
            // switching to THP is prevented if the well will
            // produce at a higher rate with THP control
            const auto& wvfpexp = this->well_ecl_.getWVFPEXP();
            bool rate_less_than_potential = true;
            if (wvfpexp.prevent()) {
                for (int p = 0; p < number_of_phases_; ++p) {
                    // Currently we use the well potentials here computed before the iterations.
                    // We may need to recompute the well potentials to get a more
                    // accurate check here.
                    rate_less_than_potential = rate_less_than_potential && (-ws.surface_rates[p]) <= ws.well_potentials[p];
                }
            }
            if(!wvfpexp.prevent() || !rate_less_than_potential) {
                this->operability_status_.thp_limit_violated_but_not_switched = false;
                return Well::ProducerCMode::THP;
            } else {
                this->operability_status_.thp_limit_violated_but_not_switched = true;
                deferred_logger.info("NOT_SWITCHING_TO_THP",
                "The THP limit is violated for producer " +
                this->name() +
                ". But the rate will increase if switched to THP. " +
                "The well is therefore kept at " + Well::ProducerCMode2String(currentControl));
            }
        }
    }

    return currentControl;
}


template <typename FluidSystem>
Well::InjectorCMode
WellInterfaceFluidSystem<FluidSystem>::
activeInjectionConstraint(const SingleWellState& ws,
                          const SummaryState& summaryState,
                          DeferredLogger& deferred_logger) const
{
    const PhaseUsage& pu = this->phaseUsage();

    const auto controls = this->well_ecl_.injectionControls(summaryState);
    const auto currentControl = ws.injection_cmode;

    if (controls.hasControl(Well::InjectorCMode::BHP) && currentControl != Well::InjectorCMode::BHP)
    {
        const auto& bhp = controls.bhp_limit;
        double current_bhp = ws.bhp;
        if (bhp < current_bhp)
            return Well::InjectorCMode::BHP;
    }

    if (controls.hasControl(Well::InjectorCMode::RATE) && currentControl != Well::InjectorCMode::RATE)
    {
        InjectorType injectorType = controls.injector_type;
        double current_rate = 0.0;

        switch (injectorType) {
        case InjectorType::WATER:
        {
            current_rate = ws.surface_rates[ pu.phase_pos[BlackoilPhases::Aqua] ];
            break;
        }
        case InjectorType::OIL:
        {
            current_rate = ws.surface_rates[ pu.phase_pos[BlackoilPhases::Liquid] ];
            break;
        }
        case InjectorType::GAS:
        {
            current_rate = ws.surface_rates[  pu.phase_pos[BlackoilPhases::Vapour] ];
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injectors " + this->well_ecl_.name());
        }

        if (controls.surface_rate < current_rate)
            return Well::InjectorCMode::RATE;
    }

    if (controls.hasControl(Well::InjectorCMode::RESV) && currentControl != Well::InjectorCMode::RESV)
    {
        double current_rate = 0.0;
        if( pu.phase_used[BlackoilPhases::Aqua] )
            current_rate += ws.reservoir_rates[ pu.phase_pos[BlackoilPhases::Aqua] ];

        if( pu.phase_used[BlackoilPhases::Liquid] )
            current_rate += ws.reservoir_rates[ pu.phase_pos[BlackoilPhases::Liquid] ];

        if( pu.phase_used[BlackoilPhases::Vapour] )
            current_rate += ws.reservoir_rates[ pu.phase_pos[BlackoilPhases::Vapour] ];

        if (controls.reservoir_rate < current_rate)
            return Well::InjectorCMode::RESV;
    }

    if (controls.hasControl(Well::InjectorCMode::THP) && currentControl != Well::InjectorCMode::THP)
    {
        const auto& thp = getTHPConstraint(summaryState);
        double current_thp = ws.thp;
        if (thp < current_thp) {
            bool rate_less_than_potential = true;
            for (int p = 0; p < number_of_phases_; ++p) {
                // Currently we use the well potentials here computed before the iterations.
                // We may need to recompute the well potentials to get a more
                // accurate check here.
                rate_less_than_potential = rate_less_than_potential && (ws.surface_rates[p]) <= ws.well_potentials[p];
            }
            if(!rate_less_than_potential) {
                this->operability_status_.thp_limit_violated_but_not_switched = false;
                return Well::InjectorCMode::THP;
            } else {
                this->operability_status_.thp_limit_violated_but_not_switched = true;
                deferred_logger.debug("NOT_SWITCHING_TO_THP",
                "The THP limit is violated for injector " +
                this->name() +
                ". But the rate will increase if switched to THP. " +
                "The well is therefore kept at " + Well::InjectorCMode2String(currentControl));
            }
        }
    }

    return currentControl;
}

template <typename FluidSystem>
bool
WellInterfaceFluidSystem<FluidSystem>::
checkIndividualConstraints(SingleWellState& ws,
                           const SummaryState& summaryState,
                           DeferredLogger& deferred_logger) const
{
    if (this->well_ecl_.isProducer()) {
        auto new_cmode = this->activeProductionConstraint(ws, summaryState, deferred_logger);
        if (new_cmode != ws.production_cmode) {
            ws.production_cmode = new_cmode;
            return true;
        }
    }

    if (this->well_ecl_.isInjector()) {
        auto new_cmode = this->activeInjectionConstraint(ws, summaryState, deferred_logger);
        if (new_cmode != ws.injection_cmode) {
            ws.injection_cmode = new_cmode;
            return true;
        }
    }

    return false;
}

template <typename FluidSystem>
std::pair<bool, double>
WellInterfaceFluidSystem<FluidSystem>::
checkGroupConstraintsInj(const Group& group,
                         const WellState& well_state,
                         const GroupState& group_state,
                         const double efficiencyFactor,
                         const Schedule& schedule,
                         const SummaryState& summaryState,
                         DeferredLogger& deferred_logger) const
{
    // Translate injector type from control to Phase.
    const auto& well_controls = this->well_ecl_.injectionControls(summaryState);
    auto injectorType = well_controls.injector_type;
    Phase injectionPhase;
    switch (injectorType) {
    case InjectorType::WATER:
    {
        injectionPhase = Phase::WATER;
        break;
    }
    case InjectorType::OIL:
    {
        injectionPhase = Phase::OIL;
        break;
    }
    case InjectorType::GAS:
    {
        injectionPhase = Phase::GAS;
        break;
    }
    default:
        throw("Expected WATER, OIL or GAS as type for injector " + name());
    }

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(phaseUsage().num_phases, 1.0);
    rateConverter_.calcInjCoeff(0, pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    const auto& ws = well_state.well(this->index_of_well_);
    // Call check for the well's injection phase.
    return WellGroupHelpers::checkGroupConstraintsInj(name(),
                                                      well_ecl_.groupName(),
                                                      group,
                                                      well_state,
                                                      group_state,
                                                      current_step_,
                                                      guide_rate_,
                                                      ws.surface_rates.data(),
                                                      injectionPhase,
                                                      phaseUsage(),
                                                      efficiencyFactor,
                                                      schedule,
                                                      summaryState,
                                                      resv_coeff,
                                                      deferred_logger);
}

template <typename FluidSystem>
std::pair<bool, double>
WellInterfaceFluidSystem<FluidSystem>::
checkGroupConstraintsProd(const Group& group,
                          const WellState& well_state,
                          const GroupState& group_state,
                          const double efficiencyFactor,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          DeferredLogger& deferred_logger) const
{
    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(this->phaseUsage().num_phases, 1.0);
    rateConverter_.calcCoeff(0, pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    const auto& ws = well_state.well(this->index_of_well_);
    return WellGroupHelpers::checkGroupConstraintsProd(name(),
                                                       well_ecl_.groupName(),
                                                       group,
                                                       well_state,
                                                       group_state,
                                                       current_step_,
                                                       guide_rate_,
                                                       ws.surface_rates.data(),
                                                       phaseUsage(),
                                                       efficiencyFactor,
                                                       schedule,
                                                       summaryState,
                                                       resv_coeff,
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
    const auto& well = well_ecl_;
    const int well_index = index_of_well_;
    auto& ws = well_state.well(well_index);

    if (well.isInjector()) {
        const auto currentControl = ws.injection_cmode;

        if (currentControl != Well::InjectorCMode::GRUP) {
            // This checks only the first encountered group limit,
            // in theory there could be several, and then we should
            // test all but the one currently applied. At that point,
            // this if-statement should be removed and we should always
            // check, skipping over only the single group parent whose
            // control is the active one for the well (if any).
            const auto& group = schedule.getGroup( well.groupName(), current_step_ );
            const double efficiencyFactor = well.getEfficiencyFactor();
            const std::pair<bool, double> group_constraint =
                checkGroupConstraintsInj(group, well_state, group_state, efficiencyFactor,
                                         schedule, summaryState, deferred_logger);
            // If a group constraint was broken, we set the current well control to
            // be GRUP.
            if (group_constraint.first) {
                ws.injection_cmode = Well::InjectorCMode::GRUP;
                const int np = well_state.numPhases();
                for (int p = 0; p<np; ++p) {
                    ws.surface_rates[p] *= group_constraint.second;
                }
            }
            return group_constraint.first;
        }
    }

    if (well.isProducer( )) {
        const auto currentControl = ws.production_cmode;

        if (currentControl != Well::ProducerCMode::GRUP) {
            // This checks only the first encountered group limit,
            // in theory there could be several, and then we should
            // test all but the one currently applied. At that point,
            // this if-statement should be removed and we should always
            // check, skipping over only the single group parent whose
            // control is the active one for the well (if any).
            const auto& group = schedule.getGroup( well.groupName(), current_step_ );
            const double efficiencyFactor = well.getEfficiencyFactor();
            const std::pair<bool, double> group_constraint =
                checkGroupConstraintsProd(group, well_state, group_state, efficiencyFactor,
                                          schedule, summaryState, deferred_logger);
            // If a group constraint was broken, we set the current well control to
            // be GRUP.
            if (group_constraint.first) {
                ws.production_cmode = Well::ProducerCMode::GRUP;
                const int np = well_state.numPhases();
                for (int p = 0; p<np; ++p) {
                    ws.surface_rates[p] *= group_constraint.second;
                }
            }
            return group_constraint.first;
        }
    }

    return false;
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
