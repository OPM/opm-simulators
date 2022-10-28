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

#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/RateConverter.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>

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
    const std::set<std::string> well_names = {"S-P2", "S-P3", "S-P4", "S-P6", "PROD1", "PROD2", "PROD3"};

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
            if (well_names.count(this->name()) > 0) {
                std::cout << " group_constraint.first (broken group control?) ? ";
                if (group_constraint.first) {
                    std::cout << " YES ";
                } else {
                    std::cout << " NO ";
                }
                std::cout << " group_constraint.second (constraint? ) : " << group_constraint.second << std::endl;
            }
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
bool
WellInterfaceFluidSystem<FluidSystem>::
checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                    const double* rates_or_potentials,
                    DeferredLogger& deferred_logger) const
{
    const PhaseUsage& pu = phaseUsage();

    if (econ_production_limits.onMinOilRate()) {
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        const double oil_rate = rates_or_potentials[pu.phase_pos[ Oil ] ];
        const double min_oil_rate = econ_production_limits.minOilRate();
        if (std::abs(oil_rate) < min_oil_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinGasRate() ) {
        assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
        const double gas_rate = rates_or_potentials[pu.phase_pos[ Gas ] ];
        const double min_gas_rate = econ_production_limits.minGasRate();
        if (std::abs(gas_rate) < min_gas_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinLiquidRate() ) {
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        const double oil_rate = rates_or_potentials[pu.phase_pos[ Oil ] ];
        const double water_rate = rates_or_potentials[pu.phase_pos[ Water ] ];
        const double liquid_rate = oil_rate + water_rate;
        const double min_liquid_rate = econ_production_limits.minLiquidRate();
        if (std::abs(liquid_rate) < min_liquid_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinReservoirFluidRate()) {
        deferred_logger.warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
    }

    return false;
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                      const SingleWellState& ws,
                      RatioLimitCheckReport& report) const
{
    assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
    assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));

    // function to calculate water cut based on rates
    auto waterCut = [](const std::vector<double>& rates,
                       const PhaseUsage& pu) {
        const double oil_rate = -rates[pu.phase_pos[Oil]];
        const double water_rate = -rates[pu.phase_pos[Water]];
        const double liquid_rate = oil_rate + water_rate;
        if (liquid_rate <= 0.)
            return 0.;
        else if (water_rate < 0)
            return 0.;
        else if (oil_rate < 0)
            return 1.;
        else
            return (water_rate / liquid_rate);

    };

    const double max_water_cut_limit = econ_production_limits.maxWaterCut();
    assert(max_water_cut_limit > 0.);

    const bool watercut_limit_violated = checkMaxRatioLimitWell(ws, max_water_cut_limit, waterCut);

    if (watercut_limit_violated) {
        report.ratio_limit_violated = true;
        checkMaxRatioLimitCompletions(ws, max_water_cut_limit, waterCut, report);
    }
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                 const SingleWellState& ws,
                 RatioLimitCheckReport& report) const
{
    assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
    assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));

    // function to calculate gor based on rates
    auto gor = [](const std::vector<double>& rates,
                  const PhaseUsage& pu) {
        const double oil_rate = -rates[pu.phase_pos[Oil]];
        const double gas_rate = -rates[pu.phase_pos[Gas]];
        if (gas_rate <= 0.)
            return 0.;
        else if (oil_rate <= 0.)
            return 1.e100; // big value to mark it as violated
        else
            return (gas_rate / oil_rate);
    };

    const double max_gor_limit = econ_production_limits.maxGasOilRatio();
    assert(max_gor_limit > 0.);

    const bool gor_limit_violated = checkMaxRatioLimitWell(ws, max_gor_limit, gor);

    if (gor_limit_violated) {
        report.ratio_limit_violated = true;
        checkMaxRatioLimitCompletions(ws, max_gor_limit, gor, report);
    }
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                 const SingleWellState& ws,
                 RatioLimitCheckReport& report) const
{
    assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
    assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));

    // function to calculate wgr based on rates
    auto wgr = [](const std::vector<double>& rates,
                  const PhaseUsage& pu) {

        const double water_rate = -rates[pu.phase_pos[Water]];
        const double gas_rate = -rates[pu.phase_pos[Gas]];
        if (water_rate <= 0.)
            return 0.;
        else if (gas_rate <= 0.)
            return 1.e100; // big value to mark it as violated
        else
            return (water_rate / gas_rate);
    };

    const double max_wgr_limit = econ_production_limits.maxWaterGasRatio();
    assert(max_wgr_limit > 0.);

    const bool wgr_limit_violated = checkMaxRatioLimitWell(ws, max_wgr_limit, wgr);

    if (wgr_limit_violated) {
        report.ratio_limit_violated = true;
        checkMaxRatioLimitCompletions(ws, max_wgr_limit, wgr, report);
    }
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                     const SingleWellState& ws,
                     RatioLimitCheckReport& report,
                     DeferredLogger& deferred_logger) const
{
    // TODO: not sure how to define the worst-offending completion when more than one
    //       ratio related limit is violated.
    //       The defintion used here is that we define the violation extent based on the
    //       ratio between the value and the corresponding limit.
    //       For each violated limit, we decide the worst-offending completion separately.
    //       Among the worst-offending completions, we use the one has the biggest violation
    //       extent.

    if (econ_production_limits.onMaxWaterCut()) {
        checkMaxWaterCutLimit(econ_production_limits, ws, report);
    }

    if (econ_production_limits.onMaxGasOilRatio()) {
        checkMaxGORLimit(econ_production_limits, ws, report);
    }

    if (econ_production_limits.onMaxWaterGasRatio()) {
        checkMaxWGRLimit(econ_production_limits, ws, report);
    }

    if (econ_production_limits.onMaxGasLiquidRatio()) {
        deferred_logger.warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
    }

    if (report.ratio_limit_violated) {
        // No worst offending completion is found because all the completions are either injecting or
        // have trivial rates.
        if(report.worst_offending_completion == INVALIDCOMPLETION) {
            std::string message = "The well ratio limit is violated but all the completion rates are trivial! " + this->name() + " is kept open";
            deferred_logger.warning("WECON_INVALIDCOMPLETION", message);
            report.ratio_limit_violated = false;
        }
        // Due to numerical instability there may exist corner cases where the well breaks
        // the ratio limit but no completion does.
        else if(report.violation_extent <= 1.) {
            std::string message = "The well ratio limit is violated but no completion ratio limit is violated! " + this->name() + " is kept open";
            deferred_logger.warning("WECON_INCONSISTANT_COMPLETION_WELL", message);
            report.ratio_limit_violated = false;
        }
    }
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
updateWellTestStateEconomic(const SingleWellState& ws,
                            const double simulation_time,
                            const bool write_message_to_opmlog,
                            WellTestState& well_test_state,
                            DeferredLogger& deferred_logger) const
{
    if (this->wellIsStopped())
        return;

    const WellEconProductionLimits& econ_production_limits = well_ecl_.getEconLimits();

    // if no limit is effective here, then continue to the next well
    if ( !econ_production_limits.onAnyEffectiveLimit() ) {
        return;
    }

    if (this->isInjector()) {
        deferred_logger.warning("ECON_LIMITS_INJECTOR_" + this->name(), this->name() + " is an injector, the production economic limits for this well will be ignored.\n");
        return;
    }

    // flag to check if the mim oil/gas rate limit is violated
    bool rate_limit_violated = false;

    const auto& quantity_limit = econ_production_limits.quantityLimit();
    if (econ_production_limits.onAnyRateLimit()) {
        if (quantity_limit == WellEconProductionLimits::QuantityLimit::POTN) {
            rate_limit_violated = checkRateEconLimits(econ_production_limits, ws.well_potentials.data(), deferred_logger);
            // Due to instability of the bhpFromThpLimit code the potentials are sometimes wrong
            // this can lead to premature shutting of wells due to rate limits of the potentials.
            // Since rates are supposed to be less or equal to the potentials, we double-check
            // that also the rate limit is violated before shutting the well.
            if (rate_limit_violated)
                rate_limit_violated = checkRateEconLimits(econ_production_limits, ws.surface_rates.data(), deferred_logger);
        }
        else {
            rate_limit_violated = checkRateEconLimits(econ_production_limits, ws.surface_rates.data(), deferred_logger);
        }
    }

    if (rate_limit_violated) {
        if (econ_production_limits.endRun()) {
            const std::string warning_message = std::string("ending run after well closed due to economic limits")
                                              + std::string("is not supported yet \n")
                                              + std::string("the program will keep running after ") + name()
                                              + std::string(" is closed");
            deferred_logger.warning("NOT_SUPPORTING_ENDRUN", warning_message);
        }

        if (econ_production_limits.validFollowonWell()) {
            deferred_logger.warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
        }

        well_test_state.close_well(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
        if (write_message_to_opmlog) {
            if (this->well_ecl_.getAutomaticShutIn()) {
                const std::string msg = std::string("well ") + name() + std::string(" will be shut due to rate economic limit");
                deferred_logger.info(msg);
            } else {
                const std::string msg = std::string("well ") + name() + std::string(" will be stopped due to rate economic limit");
                deferred_logger.info(msg);
            }
        }
        // the well is closed, not need to check other limits
        return;
    }


    if ( !econ_production_limits.onAnyRatioLimit() ) {
        // there is no need to check the ratio limits
        return;
    }

    // checking for ratio related limits, mostly all kinds of ratio.
    RatioLimitCheckReport ratio_report;

    checkRatioEconLimits(econ_production_limits, ws, ratio_report, deferred_logger);

    if (ratio_report.ratio_limit_violated) {
        const auto workover = econ_production_limits.workover();
        switch (workover) {
        case WellEconProductionLimits::EconWorkover::CON:
            {
                const int worst_offending_completion = ratio_report.worst_offending_completion;

                well_test_state.close_completion(name(), worst_offending_completion, simulation_time);
                if (write_message_to_opmlog) {
                    if (worst_offending_completion < 0) {
                        const std::string msg = std::string("Connection ") + std::to_string(- worst_offending_completion)
                                + std::string(" for well ") + name() + std::string(" will be closed due to economic limit");
                        deferred_logger.info(msg);
                    } else {
                        const std::string msg = std::string("Completion ") + std::to_string(worst_offending_completion)
                                + std::string(" for well ") + name() + std::string(" will be closed due to economic limit");
                        deferred_logger.info(msg);
                    }
                }

                bool allCompletionsClosed = true;
                const auto& connections = well_ecl_.getConnections();
                for (const auto& connection : connections) {
                    if (connection.state() == Connection::State::OPEN
                        && !well_test_state.completion_is_closed(name(), connection.complnum())) {
                        allCompletionsClosed = false;
                    }
                }

                if (allCompletionsClosed) {
                    well_test_state.close_well(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
                    if (write_message_to_opmlog) {
                        if (this->well_ecl_.getAutomaticShutIn()) {
                            const std::string msg = name() + std::string(" will be shut due to last completion closed");
                            deferred_logger.info(msg);
                        } else {
                            const std::string msg = name() + std::string(" will be stopped due to last completion closed");
                            deferred_logger.info(msg);
                        }
                    }
                }
                break;
            }
        case WellEconProductionLimits::EconWorkover::WELL:
            {
            well_test_state.close_well(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
            if (write_message_to_opmlog) {
                if (well_ecl_.getAutomaticShutIn()) {
                    // tell the control that the well is closed
                    const std::string msg = name() + std::string(" will be shut due to ratio economic limit");
                    deferred_logger.info(msg);
                } else {
                    const std::string msg = name() + std::string(" will be stopped due to ratio economic limit");
                    deferred_logger.info(msg);
                }
            }
                break;
            }
        case WellEconProductionLimits::EconWorkover::NONE:
            break;
            default:
            {
                deferred_logger.warning("NOT_SUPPORTED_WORKOVER_TYPE",
                                        "not supporting workover type " + WellEconProductionLimits::EconWorkover2String(workover) );
            }
        }
    }
}

template<typename FluidSystem>
void
WellInterfaceFluidSystem<FluidSystem>::
updateWellTestState(const SingleWellState& ws,
                    const double& simulationTime,
                    const bool& writeMessageToOPMLog,
                    WellTestState& wellTestState,
                    DeferredLogger& deferred_logger) const
{
    // updating well test state based on physical (THP/BHP) limits.
    updateWellTestStatePhysical(simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);

    // updating well test state based on Economic limits for operable wells
    if (this->isOperableAndSolvable())
        updateWellTestStateEconomic(ws, simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);

    // TODO: well can be shut/closed due to other reasons
}

template<typename FluidSystem>
template <typename RatioFunc>
void WellInterfaceFluidSystem<FluidSystem>::
checkMaxRatioLimitCompletions(const SingleWellState& ws,
                              const double max_ratio_limit,
                              const RatioFunc& ratioFunc,
                              RatioLimitCheckReport& report) const
{
    int worst_offending_completion = INVALIDCOMPLETION;

    // the maximum water cut value of the completions
    // it is used to identify the most offending completion
    double max_ratio_completion = 0;
    const int np = number_of_phases_;

    const auto& perf_data = ws.perf_data;
    const auto& perf_phase_rates = perf_data.phase_rates;
    // look for the worst_offending_completion
    for (const auto& completion : completions_) {
        std::vector<double> completion_rates(np, 0.0);

        // looping through the connections associated with the completion
        const std::vector<int>& conns = completion.second;
        for (const int c : conns) {
            for (int p = 0; p < np; ++p) {
                const double connection_rate = perf_phase_rates[c * np + p];
                completion_rates[p] += connection_rate;
            }
        } // end of for (const int c : conns)

        parallel_well_info_.communication().sum(completion_rates.data(), completion_rates.size());
        const double ratio_completion = ratioFunc(completion_rates, phaseUsage());

        if (ratio_completion > max_ratio_completion) {
            worst_offending_completion = completion.first;
            max_ratio_completion = ratio_completion;
        }
    } // end of for (const auto& completion : completions_)

    const double violation_extent = max_ratio_completion / max_ratio_limit;

    if (violation_extent > report.violation_extent) {
        report.worst_offending_completion = worst_offending_completion;
        report.violation_extent = violation_extent;
    }
}

template<typename FluidSystem>
template<typename RatioFunc>
bool WellInterfaceFluidSystem<FluidSystem>::
checkMaxRatioLimitWell(const SingleWellState& ws,
                       const double max_ratio_limit,
                       const RatioFunc& ratioFunc) const
{
    const int np = number_of_phases_;

    std::vector<double> well_rates(np, 0.0);
    for (int p = 0; p < np; ++p) {
        well_rates[p] = ws.surface_rates[p];
    }

    const double well_ratio = ratioFunc(well_rates, phaseUsage());
    return (well_ratio > max_ratio_limit);
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
    // Setting some defaults to silence warnings below.
    // Will be overwritten in the switch statement.
    Phase injectionPhase = Phase::WATER;
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
        // Should not be here.
        assert(false);
    }

    auto currentGroupControl = group_state.injection_control(group.name(), injectionPhase);
    if (currentGroupControl == Group::InjectionCMode::FLD ||
        currentGroupControl == Group::InjectionCMode::NONE) {
        if (!group.injectionGroupControlAvailable(injectionPhase)) {
            // We cannot go any further up the hierarchy. This could
            // be the FIELD group, or any group for which this has
            // been set in GCONINJE or GCONPROD. If we are here
            // anyway, it is likely that the deck set inconsistent
            // requirements, such as GRUP control mode on a well with
            // no appropriate controls defined on any of its
            // containing groups. We will therefore use the wells' bhp
            // limit equation as a fallback.
            return std::nullopt;
        } else {
            // Inject share of parents control
            const auto& parent = schedule.getGroup( group.parent(), currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            return getGroupInjectionTargetRate(parent, well_state, group_state, schedule, summaryState, injectorType, efficiencyFactor, deferred_logger);
        }
    }

    const auto pu = phaseUsage();

    if (!group.isInjectionGroup()) {
        return std::nullopt;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(pu.num_phases, 1.0);
    rateConverter_.calcCoeff(0, pvtRegionIdx(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    double sales_target = 0;
    if (schedule[currentStep()].gconsale().has(group.name())) {
        const auto& gconsale = schedule[currentStep()].gconsale().get(group.name(), summaryState);
        sales_target = gconsale.sales_target;
    }
    WellGroupHelpers::InjectionTargetCalculator tcalc(currentGroupControl, pu, resv_coeff, group.name(), sales_target, group_state, injectionPhase, group.has_gpmaint_control(injectionPhase, currentGroupControl), deferred_logger);
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, currentStep(), guideRate(), tcalc.guideTargetMode(), pu, false, injectionPhase);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child); //Note child needs to be passed to always include since the global isGrup map is not updated yet.
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.injection_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    const double orig_target = tcalc.groupTarget(group.injectionControls(injectionPhase, summaryState), deferred_logger);
    const auto chain = WellGroupHelpers::groupChainTopBot(name(), group.name(), schedule, currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || guideRate()->has(chain[ii], injectionPhase)) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    return std::max(0.0, target / efficiencyFactor);
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
    const Group::ProductionCMode& currentGroupControl = group_state.production_control(group.name());
    if (currentGroupControl == Group::ProductionCMode::FLD ||
        currentGroupControl == Group::ProductionCMode::NONE) {
        if (!group.productionGroupControlAvailable()) {
            return 1.0;
        } else {
            // Produce share of parents control
            const auto& parent = schedule.getGroup(group.parent(), currentStep());
            efficiencyFactor *= group.getGroupEfficiencyFactor();
            return getGroupProductionTargetRate(parent, well_state, group_state, schedule, summaryState, efficiencyFactor);
        }
    }

    const auto pu = phaseUsage();

    if (!group.isProductionGroup()) {
        return 1.0;
    }

    // If we are here, we are at the topmost group to be visited in the recursion.
    // This is the group containing the control we will check against.

    // Make conversion factors for RESV <-> surface rates.
    std::vector<double> resv_coeff(phaseUsage().num_phases, 1.0);
    rateConverter_.calcCoeff(0, pvtRegionIdx(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    // gconsale may adjust the grat target.
    // the adjusted rates is send to the targetCalculator
    double gratTargetFromSales = 0.0;
    if (group_state.has_grat_sales_target(group.name()))
        gratTargetFromSales = group_state.grat_sales_target(group.name());

    WellGroupHelpers::TargetCalculator tcalc(currentGroupControl, pu, resv_coeff, gratTargetFromSales, group.name(), group_state, group.has_gpmaint_control(currentGroupControl));
    WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, currentStep(), guideRate(), tcalc.guideTargetMode(), pu, true, Phase::OIL);

    auto localFraction = [&](const std::string& child) {
        return fcalc.localFraction(child, child); //Note child needs to be passed to always include since the global isGrup map is not updated yet.
    };

    auto localReduction = [&](const std::string& group_name) {
        const std::vector<double>& groupTargetReductions = group_state.production_reduction_rates(group_name);
        return tcalc.calcModeRateFromRates(groupTargetReductions);
    };

    const double orig_target = tcalc.groupTarget(group.productionControls(summaryState));
    const auto chain = WellGroupHelpers::groupChainTopBot(name(), group.name(), schedule, currentStep());
    // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
    const size_t num_ancestors = chain.size() - 1;
    double target = orig_target;
    for (size_t ii = 0; ii < num_ancestors; ++ii) {
        if ((ii == 0) || guideRate()->has(chain[ii])) {
            // Apply local reductions only at the control level
            // (top) and for levels where we have a specified
            // group guide rate.
            target -= localReduction(chain[ii]);
        }
        target *= localFraction(chain[ii+1]);
    }
    // Avoid negative target rates coming from too large local reductions.
    const double target_rate = std::max(0.0, target / efficiencyFactor);
    const auto& ws = well_state.well(this->index_of_well_);
    const auto& rates = ws.surface_rates;
    const auto current_rate = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
    double scale = 1.0;
    if (target_rate == 0.0) {
        return 0.0;
    }

    if (current_rate > 1e-14)
        scale = target_rate/current_rate;
    return scale;
}
template class WellInterfaceFluidSystem<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>>;

} // namespace Opm
