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

#include <opm/simulators/utils/DeferredLogger.hpp>

#include "opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp"
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

namespace Opm {

template<typename Scalar, typename IndexTraits>
bool WellConstraints<Scalar, IndexTraits>::
checkIndividualConstraints(SingleWellState<Scalar, IndexTraits>& ws,
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

template<typename Scalar, typename IndexTraits>
Well::InjectorCMode
WellConstraints<Scalar, IndexTraits>::
activeInjectionConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                          const SummaryState& summaryState,
                          bool& thp_limit_violated_but_not_switched,
                          DeferredLogger& deferred_logger,
                          const std::optional<Well::InjectionControls>& inj_controls) const
{
    const auto& pu = well_.phaseUsage();

    const auto controls =  inj_controls.has_value() ? inj_controls.value() : well_.wellEcl().injectionControls(summaryState);
    const auto currentControl = ws.injection_cmode;

    if (controls.hasControl(Well::InjectorCMode::BHP) && currentControl != Well::InjectorCMode::BHP)
    {
        const auto& bhp = controls.bhp_limit;
        Scalar current_bhp = ws.bhp;
        if (bhp < current_bhp)
            return Well::InjectorCMode::BHP;
    }

    if (controls.hasControl(Well::InjectorCMode::RATE) && currentControl != Well::InjectorCMode::RATE)
    {
        InjectorType injectorType = controls.injector_type;
        Scalar current_rate = 0.0;

        switch (injectorType) {
        case InjectorType::WATER:
        {
            const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            current_rate = ws.surface_rates[water_pos];
            break;
        }
        case InjectorType::OIL:
        {
            const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            current_rate = ws.surface_rates[oil_pos];
            break;
        }
        case InjectorType::GAS:
        {
            const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
            current_rate = ws.surface_rates[gas_pos];
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
        Scalar current_rate = 0.0;
        if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            current_rate += ws.reservoir_rates[water_pos];
        }

        if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
            const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            current_rate += ws.reservoir_rates[oil_pos];
        }

        if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
            const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
            current_rate += ws.reservoir_rates[gas_pos];
        }

        if (controls.reservoir_rate < current_rate)
            return Well::InjectorCMode::RESV;
    }

    // Note: we are not working on injecting network yet, so it is possible we need to change the following line
    // to be as follows to incorporate the injecting network nodal pressure
    // if (well_.wellHasTHPConstraints(summaryState) && currentControl != Well::InjectorCMode::THP)
    if (controls.hasControl(Well::InjectorCMode::THP) && currentControl != Well::InjectorCMode::THP)
    {
        const auto& thp = well_.getTHPConstraint(summaryState);
        Scalar current_thp = ws.thp;
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

template<typename Scalar, typename IndexTraits>
Well::ProducerCMode
WellConstraints<Scalar, IndexTraits>::
activeProductionConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                           const SummaryState& summaryState,
                           const RateConvFunc& calcReservoirVoidageRates,
                           bool& thp_limit_violated_but_not_switched,
                           DeferredLogger& deferred_logger,
                           const std::optional<Well::ProductionControls>& prod_controls) const
{
    const auto& pu = well_.phaseUsage();
    const auto controls = prod_controls.has_value() ? prod_controls.value() : well_.wellEcl().productionControls(summaryState);
    const auto currentControl = ws.production_cmode;

    if (controls.hasControl(Well::ProducerCMode::BHP) && currentControl != Well::ProducerCMode::BHP) {
        const Scalar bhp_limit = controls.bhp_limit;
        Scalar current_bhp = ws.bhp;
        if (bhp_limit > current_bhp)
            return Well::ProducerCMode::BHP;
    }

    if (controls.hasControl(Well::ProducerCMode::ORAT) && currentControl != Well::ProducerCMode::ORAT) {
        const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        Scalar current_rate = -ws.surface_rates[oil_pos];
        if (controls.oil_rate < current_rate)
            return Well::ProducerCMode::ORAT;
    }

    if (controls.hasControl(Well::ProducerCMode::WRAT) && currentControl != Well::ProducerCMode::WRAT) {
        const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        Scalar current_rate = -ws.surface_rates[water_pos];
        if (controls.water_rate < current_rate)
            return Well::ProducerCMode::WRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::GRAT) && currentControl != Well::ProducerCMode::GRAT) {
        const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        Scalar current_rate = -ws.surface_rates[gas_pos];
        if (controls.gas_rate < current_rate)
            return Well::ProducerCMode::GRAT;
    }

    if (controls.hasControl(Well::ProducerCMode::LRAT) && currentControl != Well::ProducerCMode::LRAT) {
        Scalar current_rate = 0.;
        if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
            const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            current_rate -= ws.surface_rates[oil_pos];
        }
        if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            current_rate -= ws.surface_rates[water_pos];
        }

        bool skip = false;
        if (controls.liquid_rate == controls.oil_rate && pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            const Scalar current_water_rate = ws.surface_rates[water_pos];
            if (std::abs(current_water_rate) < 1e-12) {
                skip = true;
                deferred_logger.debug("LRAT_ORAT_WELL", "Well " + well_.name() + " The LRAT target is equal the ORAT target and the water rate is zero, skip checking LRAT");
            }
        }
        if (!skip && controls.liquid_rate < current_rate) {
            return Well::ProducerCMode::LRAT;
        }
    }

    if (controls.hasControl(Well::ProducerCMode::RESV) && currentControl != Well::ProducerCMode::RESV) {
        Scalar current_rate = 0.0;
        if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
            const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
            current_rate -= ws.reservoir_rates[water_pos];
        }
        if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
            const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
            current_rate -= ws.reservoir_rates[oil_pos];
        }
        if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
            const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
            current_rate -= ws.reservoir_rates[gas_pos];
        }

        if (controls.prediction_mode && controls.resv_rate < current_rate)
            return Well::ProducerCMode::RESV;

        if (!controls.prediction_mode) {
            const int fipreg = 0; // not considering the region for now
            const int np = well_.numPhases();

            std::vector<Scalar> surface_rates(np, 0.0);
            if (pu.phaseIsActive(IndexTraits::waterPhaseIdx))
                surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] = controls.water_rate;
            if (pu.phaseIsActive(IndexTraits::oilPhaseIdx))
                surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] = controls.oil_rate;
            if (pu.phaseIsActive(IndexTraits::gasPhaseIdx))
                surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] = controls.gas_rate;

            std::vector<Scalar> voidage_rates(np, 0.0);
            calcReservoirVoidageRates(fipreg, well_.pvtRegionIdx(), surface_rates, voidage_rates);

            Scalar resv_rate = 0.0;
            for (int p = 0; p < np; ++p)
                resv_rate += voidage_rates[p];

            if (resv_rate < current_rate)
                return Well::ProducerCMode::RESV;
        }
    }

    if (well_.wellHasTHPConstraints(summaryState) && currentControl != Well::ProducerCMode::THP) {
        const auto& thp = well_.getTHPConstraint(summaryState);
        Scalar current_thp = ws.thp;
        // For trivial group targets (for instance caused by NETV) we dont want to flip to THP control.
        const bool dont_check = (currentControl == Well::ProducerCMode::GRUP && ws.trivial_group_target);
        if (thp > current_thp && !dont_check) {
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

template<typename Scalar, typename IndexTraits>
bool
WellConstraints<Scalar, IndexTraits>::
isFeasibleProductionConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                               const std::vector<Scalar>& surface_fractions,
                               const Well::ProducerCMode cmode_input) const
{
    auto cmode = cmode_input;
    if (cmode_input == Well::ProducerCMode::GRUP) {
        if (ws.production_cmode_group_translated.has_value()) {
            cmode = ws.production_cmode_group_translated.value();
        } else {
            return true; // can't know
        }
    }
    // We only need to check modes that can result in small fractions
    // i.e. ORAT, WRAT, GRAT, LRAT
    const bool do_check = (cmode == Well::ProducerCMode::ORAT ||
                           cmode == Well::ProducerCMode::WRAT ||
                           cmode == Well::ProducerCMode::GRAT ||
                           cmode == Well::ProducerCMode::LRAT);
    if (!do_check) {
        return true;
    }
    // somewhat ad-hoc tolerance
    Scalar fraction_tolerance = 1e-6;
    const auto& pu = well_.phaseUsage();
    // Compute fraction corresponding to the control mode
    Scalar cmode_frac = 0.0;
    if (cmode == Well::ProducerCMode::ORAT || cmode == Well::ProducerCMode::LRAT) {
        const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        cmode_frac += surface_fractions[oil_pos];
    }
    if (cmode == Well::ProducerCMode::WRAT || cmode == Well::ProducerCMode::LRAT) {
        const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        cmode_frac += surface_fractions[water_pos];
    }
    if (cmode == Well::ProducerCMode::GRAT) {
        const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        cmode_frac += surface_fractions[gas_pos];
        fraction_tolerance *= 100;
    }
    // divide by sum just in case our fractions are not normalized
    cmode_frac /= std::accumulate(surface_fractions.begin(), surface_fractions.end(), 0.0);

    return cmode_frac > fraction_tolerance;
}

template<typename Scalar, typename IndexTraits>
bool
WellConstraints<Scalar, IndexTraits>::
updateProducerControlMode(SingleWellState<Scalar, IndexTraits>& ws,
                          const SummaryState& summaryState,
                          const RateConvFunc& calcReservoirVoidageRates,
                          const Well::ProductionControls& controls,
                          const std::vector<Scalar>& surface_fractions,
                          DeferredLogger& deferred_logger) const
{
    // We first estimate strictest individual/group rate constraints, and return true if this is violated. Otherwise,
    // check if thp or bhp constraints are violated.
    // NOTE: if a well is under group control, but it's current group target is stricter than previously, this function
    // will return true (even though the mode does not change).
    const auto currentControl = ws.production_cmode;
    // use a relative (ad-hoc) tolerance
    const Scalar tol = 1e-5;

    // check most strict rate control.
    const auto [most_strict_control, most_strict_scale] =
        estimateStrictestProductionRateConstraint(ws, calcReservoirVoidageRates, controls,
                                                  surface_fractions, /*skip_zero_rate_constraints*/ false);
    if (most_strict_control != Well::ProducerCMode::CMODE_UNDEFINED && most_strict_scale < 1.0 - tol) {
        ws.production_cmode = most_strict_control;
        return true;
    }
    // check thp
    if (well_.wellHasTHPConstraints(summaryState) && currentControl != Well::ProducerCMode::THP) {
        const auto& thp = well_.getTHPConstraint(summaryState);
        Scalar current_thp = ws.thp;
        // For trivial group targets (for instance caused by NETV) we dont want to flip to THP control.
        const bool dont_check = (currentControl == Well::ProducerCMode::GRUP && ws.trivial_group_target);
        if (thp > current_thp*(1.0 + tol) && !dont_check) {
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
                //thp_limit_violated_but_not_switched = false;
                ws.production_cmode = Well::ProducerCMode::THP;
                return true;
            } else {
                //thp_limit_violated_but_not_switched = true;
                deferred_logger.info("NOT_SWITCHING_TO_THP",
                "The THP limit is violated for producer " +
                well_.name() +
                ". But the rate will increase if switched to THP. " +
                "The well is therefore kept at " + WellProducerCMode2String(currentControl));
            }
        }
    }
    // check bhp
    if (controls.hasControl(Well::ProducerCMode::BHP) && currentControl != Well::ProducerCMode::BHP) {
        const Scalar bhp_limit = controls.bhp_limit;
        Scalar current_bhp = ws.bhp;
        if (bhp_limit > current_bhp*(1.0 + tol)) {
            ws.production_cmode = Well::ProducerCMode::BHP;
            return true;
        }
    }
    // no constraints violated
    return false;
}

template<typename Scalar, typename IndexTraits>
std::pair<Well::ProducerCMode, Scalar>
WellConstraints<Scalar, IndexTraits>::
estimateStrictestProductionConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                                       const SummaryState& summaryState,
                                       const RateConvFunc& calcReservoirVoidageRates,
                                       const Well::ProductionControls& controls,
                                       const std::vector<Scalar>& surface_fractions,
                                       const bool skip_zero_rate_constraints,
                                       DeferredLogger& deferred_logger,
                                       const std::optional<Scalar> bhp_at_thp_limit) const
{
    // Estimate the most strict constraint + corresponding scaling based on current rate fractions:
    // 1. If bhp_at_thp_limit not given: potential (if available) is used to approximate pressure constraint
    //    rate-scaling - intended for initial guess
    // 2. If bhp_at_thp_limit given: we assume a converged well-state with valid ipr - intended for
    //    use within operability estimates

    const auto rates = ws.surface_rates;
    const auto tot_rates = std::accumulate(rates.begin(), rates.end(), 0.0);
    if (std::abs(tot_rates) == 0.0) {
        deferred_logger.debug("estimateStrictestProductionControl: current surface rates for well " +
                              ws.name + " are zero. Cannot determine most strict control.");
        return std::make_pair(Well::ProducerCMode::CMODE_UNDEFINED, 1.0);
    }
    Well::ProducerCMode most_strict_control = Well::ProducerCMode::CMODE_UNDEFINED;
    Scalar most_strict_scale = std::numeric_limits<Scalar>::max();

    if (!bhp_at_thp_limit.has_value() && false) {
        const auto tot_potential = std::accumulate(ws.well_potentials.begin(), ws.well_potentials.end(), 0.0);
        if (std::abs(tot_potential) > 0.0) {
            most_strict_scale = -tot_potential/tot_rates;
            if (well_.wellHasTHPConstraints(summaryState)) {
                // not neccessarily true, but most likely
                most_strict_control = Well::ProducerCMode::THP;
            } else {
                most_strict_control = Well::ProducerCMode::BHP;
            }
        }
    } else if (bhp_at_thp_limit.has_value()){
        if (*bhp_at_thp_limit > static_cast<Scalar>(controls.bhp_limit)) {
            most_strict_control = Well::ProducerCMode::THP;
        } else {
            most_strict_control = Well::ProducerCMode::BHP;
        }
        const Scalar most_strict_bhp = std::max(*bhp_at_thp_limit, static_cast<Scalar>(controls.bhp_limit));
        const Scalar tot_ipr_b = std::accumulate(ws.implicit_ipr_b.begin(), ws.implicit_ipr_b.end(), 0.0);
        const Scalar tot_ipr_a = std::accumulate(ws.implicit_ipr_a.begin(), ws.implicit_ipr_a.end(), 0.0);
        const Scalar tot_rate_at_bhp = tot_ipr_b*most_strict_bhp - tot_ipr_a;
        deferred_logger.debug("estimateStrictestProductionControl: Well " + ws.name + 
                              " has rate at bhp limit " + std::to_string(tot_rate_at_bhp) +
                              " (bhp limit " + std::to_string(most_strict_bhp) + ")" + 
                              " and current total rate " + std::to_string(tot_rates) + 
                              " ipr_a " + std::to_string(tot_ipr_a) +
                              " ipr_b " + std::to_string(tot_ipr_b));
        most_strict_scale = tot_rate_at_bhp/tot_rates;
    }
    // check rate constraints
    const auto [most_strict_rate_control, most_strict_rate_scale] =
        estimateStrictestProductionRateConstraint(ws, calcReservoirVoidageRates, controls, surface_fractions, skip_zero_rate_constraints);
    if (most_strict_rate_control != Well::ProducerCMode::CMODE_UNDEFINED && most_strict_rate_scale < most_strict_scale) {
        most_strict_scale = most_strict_rate_scale;
        most_strict_control = most_strict_rate_control;
    }
    // If we are computing e.g., well potentials, we may still have CMODE_UNDEFINED at this point. In that case, we scale according to the 
    // previous rate and keep the current control mode.
    if (most_strict_control == Well::ProducerCMode::CMODE_UNDEFINED) {
        const auto tot_prev_rates = std::accumulate(ws.prev_surface_rates.begin(), ws.prev_surface_rates.end(), 0.0);
        if (std::abs(tot_prev_rates) > 0.0) {
            most_strict_scale = std::abs(tot_prev_rates/tot_rates);
            most_strict_control = ws.production_cmode;
        } else {
            deferred_logger.debug("estimateStrictestProductionControl: previous surface rates for well " +
                                  ws.name + " are zero. This is a BUG!.");
            deferred_logger.debug("  Previous rates are " + std::to_string(ws.prev_surface_rates[0]) + ", " + std::to_string(ws.prev_surface_rates[1]) + ", " + std::to_string(ws.prev_surface_rates[2]) + ". ");
        }
    }
    return std::make_pair(most_strict_control, most_strict_scale);
}

template<typename Scalar, typename IndexTraits>
std::pair<Well::ProducerCMode, Scalar>
WellConstraints<Scalar, IndexTraits>::
estimateStrictestProductionRateConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                                          const RateConvFunc& calcReservoirVoidageRates,
                                          const Well::ProductionControls& controls,
                                          const std::vector<Scalar>& surface_fractions,
                                          const bool skip_zero_rate_constraints) const
{
    // Estimate the most strict rate constraint + corresponding scaling based on current rate fractions
    const auto rates = ws.surface_rates;
    // If all rates are zero, we can estimate strictest mode from fractions, but no scale
    const Scalar zero_rates = std::accumulate(rates.begin(), rates.end(), 0.0) == 0.0;

    Well::ProducerCMode most_strict_control = Well::ProducerCMode::CMODE_UNDEFINED;
    Scalar most_strict_scale = std::numeric_limits<Scalar>::max();

    // start with individual rate constraints
    const std::array<Well::ProducerCMode, 5> rate_modes = {Well::ProducerCMode::ORAT,
                                                          Well::ProducerCMode::WRAT,
                                                          Well::ProducerCMode::GRAT,
                                                          Well::ProducerCMode::LRAT,
                                                          Well::ProducerCMode::RESV};
    for (const auto& mode : rate_modes) {
        if (!controls.hasControl(mode))
            continue;
        const Scalar scale = getProductionControlModeScale(ws, calcReservoirVoidageRates, mode, controls, surface_fractions, skip_zero_rate_constraints);
        if (scale >= 0.0 && scale < most_strict_scale) {
            most_strict_scale = scale;
            most_strict_control = mode;
        }
    }
    // check group constraints if target is given in well-state
    if (controls.hasControl(Well::ProducerCMode::GRUP) && !ws.prevent_group_control && ws.group_target.has_value() && ws.production_cmode_group_translated.has_value()) {
        const Scalar scale = getProductionControlModeScale(ws, calcReservoirVoidageRates, ws.production_cmode_group_translated.value(), controls,
                                                           surface_fractions, skip_zero_rate_constraints, ws.group_target.value());
        if (scale >= 0.0 && scale < most_strict_scale) {
            most_strict_scale = scale;
            most_strict_control = Well::ProducerCMode::GRUP;
        }
    }
    if (zero_rates) {
        most_strict_scale = std::numeric_limits<Scalar>::max();
    }
    return std::make_pair(most_strict_control, most_strict_scale);
}

template<typename Scalar, typename IndexTraits>
Scalar
WellConstraints<Scalar, IndexTraits>::
getProductionControlModeScale(const SingleWellState<Scalar, IndexTraits>& ws,
                              const RateConvFunc& calcReservoirVoidageRates,
                              const Well::ProducerCMode& cmode,
                              const Well::ProductionControls& controls,
                              const std::vector<Scalar>& surface_fractions,
                              const bool skip_zero_rate_constraints,
                              const std::optional<Scalar> target) const
{
    // check if cmode is feasible with current rate fractions
    if (!isFeasibleProductionConstraint(ws, surface_fractions, cmode)) {
        return -1.0;
    }
    const auto& pu = well_.phaseUsage();
    Scalar current_rate = 0.0;
    Scalar target_rate = 0.0;
    // use fractions if all rates are zero 
    const bool zero_rates = std::accumulate(ws.surface_rates.begin(), ws.surface_rates.end(), 0.0) == 0.0;
    const auto rates = zero_rates ? surface_fractions : ws.surface_rates;
    switch (cmode) {
        case Well::ProducerCMode::ORAT:
            current_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
            target_rate = target.has_value() ? *target : controls.oil_rate;
            break;
        case Well::ProducerCMode::WRAT:
            current_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
            target_rate = target.has_value() ? *target : controls.water_rate;
            break;
        case Well::ProducerCMode::GRAT:
            current_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)];
            target_rate = target.has_value() ? *target : controls.gas_rate;
            break;
        case Well::ProducerCMode::LRAT:
            current_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
            current_rate += -rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
            target_rate = target.has_value() ? *target : controls.liquid_rate;
            break;
        case Well::ProducerCMode::RESV:
            for (int p = 0; p < well_.numPhases(); ++p) {
                current_rate += -ws.reservoir_rates[p];
            }
            if (controls.prediction_mode || target.has_value()) {
                target_rate = target.has_value() ? *target : controls.resv_rate;
            } else {
                const int fipreg = 0; // not considering the region for now
                const int np = well_.numPhases();

                std::vector<Scalar> target_surface_rates(np, 0.0);
                if (pu.phaseIsActive(IndexTraits::waterPhaseIdx))
                    target_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] = controls.water_rate;
                if (pu.phaseIsActive(IndexTraits::oilPhaseIdx))
                    target_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] = controls.oil_rate;
                if (pu.phaseIsActive(IndexTraits::gasPhaseIdx))
                    target_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] = controls.gas_rate;

                std::vector<Scalar> target_voidage_rates(np, 0.0);
                calcReservoirVoidageRates(fipreg, well_.pvtRegionIdx(), target_surface_rates, target_voidage_rates);
                // Scalar resv_rate = 0.0;
                for (int p = 0; p < np; ++p) {
                    target_rate += target_voidage_rates[p];
                }
            }
            break;
        default:
            // undefined mode, no scaling applied
            return -1.0;
            break;
    }
    const bool valid_scale = (!skip_zero_rate_constraints || std::abs(target_rate) > 0.0);
    return  valid_scale ? std::abs(target_rate/current_rate) : -1.0;
}

template class WellConstraints<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class WellConstraints<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
