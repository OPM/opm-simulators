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

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>

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
            OPM_THROW(std::runtime_error,
                      "Expected WATER, OIL or GAS as type for injectors " + well_.name());
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
Scalar
WellConstraints<Scalar, IndexTraits>::
getProductionControlModeScale(const std::vector<Scalar>& pos_surface_rates,
                              const std::vector<Scalar>& pos_reservoir_rates,
                              const RateConvFunc& calcReservoirVoidageRates,
                              const Well::ProducerCMode& cmode,
                              const Well::ProductionControls& controls,
                              const std::optional<Scalar> target) const
{
    // pos_surface_rates[p]  > 0 for production (sign-flipped from ws.surface_rates)
    // pos_reservoir_rates[p] > 0 for production (sign-flipped from ws.reservoir_rates)
    const auto& pu = well_.phaseUsage();
    Scalar current_rate = Scalar(0);
    Scalar target_rate  = Scalar(0);
    switch (cmode) {
        case Well::ProducerCMode::ORAT:
            current_rate = pos_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
            target_rate  = target.value_or(controls.oil_rate);
            break;
        case Well::ProducerCMode::WRAT:
            current_rate = pos_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
            target_rate  = target.value_or(controls.water_rate);
            break;
        case Well::ProducerCMode::GRAT:
            current_rate = pos_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)];
            target_rate  = target.value_or(controls.gas_rate);
            break;
        case Well::ProducerCMode::LRAT:
            current_rate = pos_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
            current_rate += pos_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
            target_rate  = target.value_or(controls.liquid_rate);
            break;
        case Well::ProducerCMode::RESV:
            for (int p = 0; p < well_.numPhases(); ++p) {
                current_rate += pos_reservoir_rates[p];
            }
            if (!controls.prediction_mode || target.has_value()) {
                target_rate = target.value_or(controls.resv_rate);
            } else {
                const int fipreg = 0;
                const int np = well_.numPhases();
                std::vector<Scalar> target_surface_rates(np, Scalar(0));
                if (pu.phaseIsActive(IndexTraits::waterPhaseIdx))
                    target_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] = controls.water_rate;
                if (pu.phaseIsActive(IndexTraits::oilPhaseIdx))
                    target_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] = controls.oil_rate;
                if (pu.phaseIsActive(IndexTraits::gasPhaseIdx))
                    target_surface_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] = controls.gas_rate;
                std::vector<Scalar> target_voidage_rates(np, Scalar(0));
                calcReservoirVoidageRates(fipreg, well_.pvtRegionIdx(), target_surface_rates, target_voidage_rates);
                for (int p = 0; p < np; ++p)
                    target_rate += target_voidage_rates[p];
            }
            break;
        default:
            break;
    }
    // If zero target, scale is always zero
    if (target_rate == Scalar(0)) return Scalar(0);

    return std::abs(current_rate) > Scalar(0) 
        ? std::abs(target_rate / current_rate) 
        : std::numeric_limits<Scalar>::max();
}

template<typename Scalar, typename IndexTraits>
std::pair<Well::ProducerCMode, Scalar>
WellConstraints<Scalar, IndexTraits>::
estimateStrictestRateConstraintFromRatesImpl_(const std::vector<Scalar>& pos_surface_rates,
                                              const std::vector<Scalar>& pos_reservoir_rates,
                                              const RateConvFunc& calcReservoirVoidageRates,
                                              const Well::ProductionControls& controls) const
{
    Well::ProducerCMode most_strict_control = Well::ProducerCMode::CMODE_UNDEFINED;
    Scalar most_strict_scale = std::numeric_limits<Scalar>::max();

    constexpr std::array<Well::ProducerCMode, 5> rate_modes = {Well::ProducerCMode::ORAT,
                                                               Well::ProducerCMode::WRAT,
                                                               Well::ProducerCMode::GRAT,
                                                               Well::ProducerCMode::LRAT,
                                                               Well::ProducerCMode::RESV};
    for (const auto& mode : rate_modes) {
        if (!controls.hasControl(mode))
            continue;
        const Scalar scale = getProductionControlModeScale(
            pos_surface_rates, pos_reservoir_rates, calcReservoirVoidageRates, mode, controls);
        if (scale >= Scalar(0) && scale < most_strict_scale) {
            most_strict_scale = scale;
            most_strict_control = mode;
        }
    }
    return {most_strict_control, most_strict_scale};
}

template<typename Scalar, typename IndexTraits>
std::pair<Well::ProducerCMode, Scalar>
WellConstraints<Scalar, IndexTraits>::
estimateStrictestProductionConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                                      const RateConvFunc& calcReservoirVoidageRates,
                                      const Well::ProductionControls& controls,
                                      const bool check_group_constraints,
                                      DeferredLogger& deferred_logger,
                                      const std::optional<Scalar> bhp_at_thp_limit) const
{
    // Returns (control_mode, scale) where scale applied to current rates gives the constrained rates.
    // bhp_at_thp_limit must be provided (from a converged IPR) for the pressure-constraint term.
    const auto& surf = ws.surface_rates;
    const auto tot_rates = std::accumulate(surf.begin(), surf.end(), Scalar(0));

    Well::ProducerCMode most_strict_control = Well::ProducerCMode::CMODE_UNDEFINED;
    Scalar most_strict_scale = std::numeric_limits<Scalar>::max();

    // Pressure constraint: requires a converged IPR
    if (bhp_at_thp_limit.has_value() && std::abs(tot_rates) > Scalar(0)) {
        most_strict_control = (*bhp_at_thp_limit > controls.bhp_limit)
                              ? Well::ProducerCMode::THP
                              : Well::ProducerCMode::BHP;
        const Scalar most_strict_bhp = std::max(*bhp_at_thp_limit, controls.bhp_limit);
        const Scalar tot_ipr_b = std::accumulate(ws.implicit_ipr_b.begin(), ws.implicit_ipr_b.end(), Scalar(0));
        const Scalar tot_ipr_a = std::accumulate(ws.implicit_ipr_a.begin(), ws.implicit_ipr_a.end(), Scalar(0));
        const Scalar tot_rate_at_bhp = tot_ipr_b * most_strict_bhp - tot_ipr_a;
        most_strict_scale = tot_rate_at_bhp / tot_rates;
    }

    // Rate constraints
    const auto [most_strict_rate_control, most_strict_rate_scale] =
        estimateStrictestProductionRateConstraint(ws, calcReservoirVoidageRates,
                                                  controls, check_group_constraints, deferred_logger);
    if (most_strict_rate_control != Well::ProducerCMode::CMODE_UNDEFINED &&
        most_strict_rate_scale < most_strict_scale) {
        most_strict_scale   = most_strict_rate_scale;
        most_strict_control = most_strict_rate_control;
    }

    // If still undefined (e.g. computing well potentials with no active pressure constraint),
    // fall back to previous-timestep rate to keep the scale finite.
    if (most_strict_control == Well::ProducerCMode::CMODE_UNDEFINED) {
        const auto tot_prev = std::accumulate(ws.prev_surface_rates.begin(), ws.prev_surface_rates.end(), Scalar(0));
        if (std::abs(tot_prev) > Scalar(0)) {
            most_strict_scale   = std::abs(tot_prev / tot_rates);
            most_strict_control = ws.production_cmode;
        } else {
            deferred_logger.debug("estimateStrictestProductionControl: previous surface rates for well " +
                                  ws.name + " are zero. This is a BUG!.");
            deferred_logger.debug("  Previous rates are " + std::to_string(ws.prev_surface_rates[0]) + ", " +
                                  std::to_string(ws.prev_surface_rates[1]) + ", " +
                                  std::to_string(ws.prev_surface_rates[2]) + ". ");
        }
    }
    return {most_strict_control, most_strict_scale};
}

template<typename Scalar, typename IndexTraits>
std::pair<Well::ProducerCMode, Scalar>
WellConstraints<Scalar, IndexTraits>::
estimateStrictestProductionRateConstraint(const SingleWellState<Scalar, IndexTraits>& ws,
                                          const RateConvFunc& calcReservoirVoidageRates,
                                          const Well::ProductionControls& controls,
                                          const bool check_group_constraints,
                                          DeferredLogger& deferred_logger) const
{
    auto surf = ws.surface_rates;
    const auto tot_rates = std::accumulate(surf.begin(), surf.end(), Scalar(0));
    if (std::abs(tot_rates) == Scalar(0)) {
        // we likely have a zero rate-constraint, use previous rates to determine the most strict control
        surf = ws.prev_surface_rates;
        const auto tot_prev_rates = std::accumulate(surf.begin(), surf.end(), Scalar(0));
        if (std::abs(tot_prev_rates) == Scalar(0)) {
            deferred_logger.debug("estimateStrictestProductionRateControl: current and previous surface rates for well " +
                                ws.name + " are zero. Cannot determine most strict control.");
        return {Well::ProducerCMode::CMODE_UNDEFINED, Scalar(1)};
        }
    }

    // Build positive-valued rate vectors (producers have negative surface_rates / reservoir_rates)
    const int np = static_cast<int>(surf.size());
    std::vector<Scalar> pos_surf(np), pos_resv(np);
    for (int p = 0; p < np; ++p) {
        pos_surf[p] = -surf[p];
        // 
        pos_resv[p] = -ws.reservoir_rates[p];  
    }

    auto [most_strict_control, most_strict_scale] =
        estimateStrictestRateConstraintFromRatesImpl_(
            pos_surf, pos_resv, calcReservoirVoidageRates, controls);

    // Group constraint check — only meaningful with actual well-state rates
    if (check_group_constraints && controls.hasControl(Well::ProducerCMode::GRUP) &&
        ws.group_target.has_value()) {
        const bool use_fallback = ws.group_target_fallback.has_value() && ws.use_group_target_fallback;
        const Scalar target = use_fallback
            ? ws.group_target_fallback.value().target_value
            : ws.group_target.value().target_value;
        const Group::ProductionCMode gmode = use_fallback
            ? ws.group_target_fallback.value().production_cmode
            : ws.group_target.value().production_cmode;
        const auto cmode = WellProducerCModeFromString(Group::ProductionCMode2String(gmode));
        const Scalar scale = getProductionControlModeScale(
            pos_surf, pos_resv, calcReservoirVoidageRates, cmode, controls, target);
        if (scale >= Scalar(0) && scale < most_strict_scale) {
            most_strict_scale   = scale;
            most_strict_control = Well::ProducerCMode::GRUP;
        }
    }
    return {most_strict_control, most_strict_scale};
}

template<typename Scalar, typename IndexTraits>
std::pair<Well::ProducerCMode, Scalar>
WellConstraints<Scalar, IndexTraits>::
estimateStrictestProductionRateConstraintFromRates(
    const std::vector<Scalar>& pos_surface_rates,
    const std::vector<Scalar>& pos_reservoir_rates,
    const RateConvFunc& calcReservoirVoidageRates,
    const Well::ProductionControls& controls,
    DeferredLogger& deferred_logger) const
{
    // Find the strictest rate constraint relative to explicit (positive-valued) rate vectors.
    // Intended for the well-potentials path where pos_surface_rates = ws.well_potentials and
    // pos_reservoir_rates = reservoir equivalents computed from the potentials.
    // GRUP constraints are not checked because potentials do not represent group allocation.
    const auto tot_rate = std::accumulate(
        pos_surface_rates.begin(), pos_surface_rates.end(), Scalar(0));
    if (tot_rate <= Scalar(0)) {
        deferred_logger.debug("estimateStrictestProductionRateConstraintFromRates: potentials for well " +
                              well_.name() + " are zero or negative.");
        return {Well::ProducerCMode::CMODE_UNDEFINED, Scalar(1)};
    }
    return estimateStrictestRateConstraintFromRatesImpl_(
        pos_surface_rates, pos_reservoir_rates, calcReservoirVoidageRates, controls);
}

template class WellConstraints<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class WellConstraints<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
