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
#include <opm/simulators/wells/WellTest.hpp>

#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

namespace Opm
{

template<class RatioFunc>
bool WellTest::checkMaxRatioLimitWell(const SingleWellState& ws,
                                      const double max_ratio_limit,
                                      const RatioFunc& ratioFunc) const
{
    const int np = well_.numPhases();

    std::vector<double> well_rates(np, 0.0);
    for (int p = 0; p < np; ++p) {
        well_rates[p] = ws.surface_rates[p];
    }

    const double well_ratio = ratioFunc(well_rates, well_.phaseUsage());
    return (well_ratio > max_ratio_limit);
}

template<class RatioFunc>
void WellTest::checkMaxRatioLimitCompletions(const SingleWellState& ws,
                                             const double max_ratio_limit,
                                             const RatioFunc& ratioFunc,
                                             RatioLimitCheckReport& report) const
{
    int worst_offending_completion = RatioLimitCheckReport::INVALIDCOMPLETION;

    // the maximum water cut value of the completions
    // it is used to identify the most offending completion
    double max_ratio_completion = 0;
    const int np = well_.numPhases();

    const auto& perf_data = ws.perf_data;
    const auto& perf_phase_rates = perf_data.phase_rates;
    // look for the worst_offending_completion
    for (const auto& completion : well_.getCompletions()) {
        std::vector<double> completion_rates(np, 0.0);

        // looping through the connections associated with the completion
        const std::vector<int>& conns = completion.second;
        for (const int c : conns) {
            for (int p = 0; p < np; ++p) {
                const double connection_rate = perf_phase_rates[c * np + p];
                completion_rates[p] += connection_rate;
            }
        } // end of for (const int c : conns)

        well_.parallelWellInfo().communication().sum(completion_rates.data(), completion_rates.size());
        const double ratio_completion = ratioFunc(completion_rates, well_.phaseUsage());

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

void WellTest::checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                                const SingleWellState& ws,
                                RatioLimitCheckReport& report) const
{
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

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

    const bool gor_limit_violated = this->checkMaxRatioLimitWell(ws, max_gor_limit, gor);

    if (gor_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_gor_limit, gor, report);
    }
}

void WellTest::checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                                const SingleWellState& ws,
                                RatioLimitCheckReport& report) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Water = BlackoilPhases::Aqua;

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

    const bool wgr_limit_violated = this->checkMaxRatioLimitWell(ws, max_wgr_limit, wgr);

    if (wgr_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_wgr_limit, wgr, report);
    }
}

void WellTest::checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                                     const SingleWellState& ws,
                                     RatioLimitCheckReport& report) const
{
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

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

    const bool watercut_limit_violated =
        this->checkMaxRatioLimitWell(ws, max_water_cut_limit, waterCut);

    if (watercut_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_water_cut_limit,
                                            waterCut, report);
    }
}

bool WellTest::checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                   const std::vector<double>& rates_or_potentials,
                                   DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    const PhaseUsage& pu = well_.phaseUsage();

    if (econ_production_limits.onMinOilRate()) {
        const double oil_rate = rates_or_potentials[pu.phase_pos[Oil]];
        const double min_oil_rate = econ_production_limits.minOilRate();
        if (std::abs(oil_rate) < min_oil_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinGasRate() ) {
        const double gas_rate = rates_or_potentials[pu.phase_pos[Gas]];
        const double min_gas_rate = econ_production_limits.minGasRate();
        if (std::abs(gas_rate) < min_gas_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinLiquidRate() ) {
        const double oil_rate = rates_or_potentials[pu.phase_pos[Oil]];
        const double water_rate = rates_or_potentials[pu.phase_pos[Water]];
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

WellTest::RatioLimitCheckReport WellTest::
checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                     const SingleWellState& ws,
                     DeferredLogger& deferred_logger) const
{
    // TODO: not sure how to define the worst-offending completion when more than one
    //       ratio related limit is violated.
    //       The defintion used here is that we define the violation extent based on the
    //       ratio between the value and the corresponding limit.
    //       For each violated limit, we decide the worst-offending completion separately.
    //       Among the worst-offending completions, we use the one has the biggest violation
    //       extent.
    RatioLimitCheckReport report;

    if (econ_production_limits.onMaxWaterCut()) {
        this->checkMaxWaterCutLimit(econ_production_limits, ws, report);
    }

    if (econ_production_limits.onMaxGasOilRatio()) {
        this->checkMaxGORLimit(econ_production_limits, ws, report);
    }

    if (econ_production_limits.onMaxWaterGasRatio()) {
        this->checkMaxWGRLimit(econ_production_limits, ws, report);
    }

    if (econ_production_limits.onMaxGasLiquidRatio()) {
        deferred_logger.warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
    }

    if (report.ratio_limit_violated) {
        // No worst offending completion is found because all the completions are either injecting or
        // have trivial rates.
        if(report.worst_offending_completion == RatioLimitCheckReport::INVALIDCOMPLETION) {
            std::string message = "The well ratio limit is violated but all the completion rates are trivial! " + well_.name() + " is kept open";
            deferred_logger.warning("WECON_INVALIDCOMPLETION", message);
            report.ratio_limit_violated = false;
        }
        // Due to numerical instability there may exist corner cases where the well breaks
        // the ratio limit but no completion does.
        else if(report.violation_extent <= 1.) {
            std::string message = "The well ratio limit is violated but no completion ratio limit is violated! " + well_.name() + " is kept open";
            deferred_logger.warning("WECON_INCONSISTANT_COMPLETION_WELL", message);
            report.ratio_limit_violated = false;
        }
    }

    return report;
}

void WellTest::updateWellTestStateEconomic(const SingleWellState& ws,
                                           const double simulation_time,
                                           const bool write_message_to_opmlog,
                                           WellTestState& well_test_state,
                                           DeferredLogger& deferred_logger) const
{
    if (well_.wellIsStopped())
        return;

    const WellEconProductionLimits& econ_production_limits = well_.wellEcl().getEconLimits();

    // if no limit is effective here, then continue to the next well
    if (!econ_production_limits.onAnyEffectiveLimit()) {
        return;
    }

    if (well_.isInjector()) {
        deferred_logger.warning("ECON_LIMITS_INJECTOR_" + well_.name(), well_.name() + " is an injector, the production economic limits for this well will be ignored.\n");
        return;
    }

    // flag to check if the mim oil/gas rate limit is violated
    bool rate_limit_violated = false;

    const auto& quantity_limit = econ_production_limits.quantityLimit();
    if (econ_production_limits.onAnyRateLimit()) {
        if (quantity_limit == WellEconProductionLimits::QuantityLimit::POTN) {
            rate_limit_violated = this->checkRateEconLimits(econ_production_limits,
                                                            ws.well_potentials,
                                                            deferred_logger);
            // Due to instability of the bhpFromThpLimit code the potentials are sometimes wrong
            // this can lead to premature shutting of wells due to rate limits of the potentials.
            // Since rates are supposed to be less or equal to the potentials, we double-check
            // that also the rate limit is violated before shutting the well.
            if (rate_limit_violated)
                rate_limit_violated = this->checkRateEconLimits(econ_production_limits,
                                                                ws.surface_rates,
                                                                deferred_logger);
        }
        else {
            rate_limit_violated = this->checkRateEconLimits(econ_production_limits,
                                                            ws.surface_rates,
                                                            deferred_logger);
        }
    }

    if (rate_limit_violated) {
        if (econ_production_limits.endRun()) {
            const std::string warning_message = std::string("ending run after well closed due to economic limits")
                                              + std::string("is not supported yet \n")
                                              + std::string("the program will keep running after ") + well_.name()
                                              + std::string(" is closed");
            deferred_logger.warning("NOT_SUPPORTING_ENDRUN", warning_message);
        }

        if (econ_production_limits.validFollowonWell()) {
            deferred_logger.warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
        }

        well_test_state.close_well(well_.name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
        if (write_message_to_opmlog) {
            if (well_.wellEcl().getAutomaticShutIn()) {
                const std::string msg = std::string("well ") + well_.name() + std::string(" will be shut due to rate economic limit");
                deferred_logger.info(msg);
            } else {
                const std::string msg = std::string("well ") + well_.name() + std::string(" will be stopped due to rate economic limit");
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
    RatioLimitCheckReport ratio_report =
        this->checkRatioEconLimits(econ_production_limits, ws, deferred_logger);

    if (ratio_report.ratio_limit_violated) {
        const auto workover = econ_production_limits.workover();
        switch (workover) {
        case WellEconProductionLimits::EconWorkover::CON:
            {
                const int worst_offending_completion = ratio_report.worst_offending_completion;

                well_test_state.close_completion(well_.name(), worst_offending_completion, simulation_time);
                if (write_message_to_opmlog) {
                    if (worst_offending_completion < 0) {
                        const std::string msg = std::string("Connection ") + std::to_string(- worst_offending_completion)
                                + std::string(" for well ") + well_.name() + std::string(" will be closed due to economic limit");
                        deferred_logger.info(msg);
                    } else {
                        const std::string msg = std::string("Completion ") + std::to_string(worst_offending_completion)
                                + std::string(" for well ") + well_.name() + std::string(" will be closed due to economic limit");
                        deferred_logger.info(msg);
                    }
                }

                bool allCompletionsClosed = true;
                const auto& connections = well_.wellEcl().getConnections();
                for (const auto& connection : connections) {
                    if (connection.state() == Connection::State::OPEN
                        && !well_test_state.completion_is_closed(well_.name(), connection.complnum())) {
                        allCompletionsClosed = false;
                    }
                }

                if (allCompletionsClosed) {
                    well_test_state.close_well(well_.name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
                    if (write_message_to_opmlog) {
                        if (well_.wellEcl().getAutomaticShutIn()) {
                            const std::string msg = well_.name() + std::string(" will be shut due to last completion closed");
                            deferred_logger.info(msg);
                        } else {
                            const std::string msg = well_.name() + std::string(" will be stopped due to last completion closed");
                            deferred_logger.info(msg);
                        }
                    }
                }
                break;
            }
        case WellEconProductionLimits::EconWorkover::WELL:
            {
            well_test_state.close_well(well_.name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
            if (write_message_to_opmlog) {
                if (well_.wellEcl().getAutomaticShutIn()) {
                    // tell the control that the well is closed
                    const std::string msg = well_.name() + std::string(" will be shut due to ratio economic limit");
                    deferred_logger.info(msg);
                } else {
                    const std::string msg = well_.name() + std::string(" will be stopped due to ratio economic limit");
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
                                        "not supporting workover type " + WellEconProductionLimits::EconWorkover2String(workover));
            }
        }
    }
}

void WellTest::updateWellTestStatePhysical(const double simulation_time,
                                           const bool write_message_to_opmlog,
                                           WellTestState& well_test_state,
                                           DeferredLogger& deferred_logger) const
{
    if (well_test_state.well_is_closed(well_.name())) {
        // Already closed, do nothing.
    } else {
        well_test_state.close_well(well_.name(), WellTestConfig::Reason::PHYSICAL, simulation_time);
        if (write_message_to_opmlog) {
            const std::string action = well_.wellEcl().getAutomaticShutIn() ? "shut" : "stopped";
            const std::string msg = "Well " + well_.name()
                + " will be " + action + " as it can not operate under current reservoir conditions.";
            deferred_logger.info(msg);
        }
    }
}

} // namespace Opm
