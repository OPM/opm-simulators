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
#include <opm/simulators/wells/WellConstraints.hpp>

#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

namespace Opm
{

bool WellConstraints::
checkIndividualConstraints(SingleWellState& ws,
                           const SummaryState& summaryState,
                           const RateConvFunc& calcReservoirVoidageRates,
                           bool& thp_limit_violated_but_not_switched,
                           DeferredLogger& deferred_logger,
                           const std::optional<Well::InjectionControls>& inj_controls,
                           const std::optional<Well::ProductionControls>& prod_controls) const
{
    if (well_.isProducer()) {
        auto new_cmode = this->activeProductionConstraint(ws, summaryState,
                                                          calcReservoirVoidageRates,
                                                          thp_limit_violated_but_not_switched,
                                                          deferred_logger,
                                                          prod_controls);
        if (new_cmode != ws.production_cmode) {
            ws.production_cmode = new_cmode;
            return true;
        }
    }

    if (well_.isInjector()) {
        auto new_cmode = this->activeInjectionConstraint(ws, summaryState,
                                                        thp_limit_violated_but_not_switched,
                                                        deferred_logger,
                                                        inj_controls);
        if (new_cmode != ws.injection_cmode) {
            ws.injection_cmode = new_cmode;
            return true;
        }
    }

    return false;
}

Well::InjectorCMode WellConstraints::
activeInjectionConstraint(const SingleWellState& ws,
                          const SummaryState& summaryState,
                          bool& thp_limit_violated_but_not_switched,
                          DeferredLogger& deferred_logger,
                          const std::optional<Well::InjectionControls>& inj_controls) const
{
    const PhaseUsage& pu = well_.phaseUsage();

    const auto controls =  inj_controls.has_value() ? inj_controls.value() : well_.wellEcl().injectionControls(summaryState);
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
            throw("Expected WATER, OIL or GAS as type for injectors " + well_.name());
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

    // Note: we are not working on injecting network yet, so it is possible we need to change the following line
    // to be as follows to incorporate the injecting network nodal pressure
    // if (well_.wellHasTHPConstraints(summaryState) && currentControl != Well::InjectorCMode::THP)
    if (controls.hasControl(Well::InjectorCMode::THP) && currentControl != Well::InjectorCMode::THP)
    {
        const auto& thp = well_.getTHPConstraint(summaryState);
        double current_thp = ws.thp;
        if (thp < current_thp) {
            bool rate_less_than_potential = true;
            for (int p = 0; p < well_.numPhases(); ++p) {
                // Currently we use the well potentials here computed before the iterations.
                // We may need to recompute the well potentials to get a more
                // accurate check here.
                rate_less_than_potential = rate_less_than_potential && (ws.surface_rates[p]) <= ws.well_potentials[p];
            }
            if (!rate_less_than_potential) {
                thp_limit_violated_but_not_switched = false;
                return Well::InjectorCMode::THP;
            } else {
                thp_limit_violated_but_not_switched = true;
                deferred_logger.debug("NOT_SWITCHING_TO_THP",
                "The THP limit is violated for injector " +
                well_.name() +
                ". But the rate will increase if switched to THP. " +
                "The well is therefore kept at " + WellInjectorCMode2String(currentControl));
            }
        }
    }

    return currentControl;
}

Well::ProducerCMode WellConstraints::
activeProductionConstraint(const SingleWellState& ws,
                           const SummaryState& summaryState,
                           const RateConvFunc& calcReservoirVoidageRates,
                           bool& thp_limit_violated_but_not_switched,
                           DeferredLogger& deferred_logger,
                           const std::optional<Well::ProductionControls>& prod_controls) const
{
    const PhaseUsage& pu = well_.phaseUsage();
    const auto controls = prod_controls.has_value() ? prod_controls.value() : well_.wellEcl().productionControls(summaryState);
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
                deferred_logger.debug("LRAT_ORAT_WELL", "Well " + well_.name() + " The LRAT target is equal the ORAT target and the water rate is zero, skip checking LRAT");
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
            const int np = well_.numPhases();

            std::vector<double> surface_rates(np, 0.0);
            if (pu.phase_used[BlackoilPhases::Aqua])
                surface_rates[pu.phase_pos[BlackoilPhases::Aqua]] = controls.water_rate;
            if (pu.phase_used[BlackoilPhases::Liquid])
                surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = controls.oil_rate;
            if (pu.phase_used[BlackoilPhases::Vapour])
                surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = controls.gas_rate;

            std::vector<double> voidage_rates(np, 0.0);
            calcReservoirVoidageRates(fipreg, well_.pvtRegionIdx(), surface_rates, voidage_rates);

            double resv_rate = 0.0;
            for (int p = 0; p < np; ++p)
                resv_rate += voidage_rates[p];

            if (resv_rate < current_rate)
                return Well::ProducerCMode::RESV;
        }
    }

    if (well_.wellHasTHPConstraints(summaryState) && currentControl != Well::ProducerCMode::THP) {
        const auto& thp = well_.getTHPConstraint(summaryState);
        double current_thp = ws.thp;
        if (thp > current_thp && !ws.trivial_target) {
            // If WVFPEXP item 4 is set to YES1 or YES2
            // switching to THP is prevented if the well will
            // produce at a higher rate with THP control
            const auto& wvfpexp = well_.wellEcl().getWVFPEXP();
            bool rate_less_than_potential = true;
            if (wvfpexp.prevent()) {
                for (int p = 0; p < well_.numPhases(); ++p) {
                    // Currently we use the well potentials here computed before the iterations.
                    // We may need to recompute the well potentials to get a more
                    // accurate check here.
                    rate_less_than_potential = rate_less_than_potential && (-ws.surface_rates[p]) <= ws.well_potentials[p];
                }
            }
            if (!wvfpexp.prevent() || !rate_less_than_potential) {
                thp_limit_violated_but_not_switched = false;
                return Well::ProducerCMode::THP;
            } else {
                thp_limit_violated_but_not_switched = true;
                deferred_logger.info("NOT_SWITCHING_TO_THP",
                "The THP limit is violated for producer " +
                well_.name() +
                ". But the rate will increase if switched to THP. " +
                "The well is therefore kept at " + WellProducerCMode2String(currentControl));
            }
        }
    }

    return currentControl;
}

} // namespace Opm
