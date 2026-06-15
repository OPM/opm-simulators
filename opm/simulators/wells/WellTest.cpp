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

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/common/utility/TimeService.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_set>

namespace {

//! \brief Calendar date (DD-Mon-YYYY) reached at start_time + sim_time seconds.
std::string dateString(const std::time_t start_time, const double sim_time)
{
    const std::time_t cur_time = Opm::TimeService::advance(start_time, sim_time);
    return fmt::format("{:%d-%b-%Y}", fmt::gmtime(cur_time));
}

std::string ceconDateString(const std::time_t start_time, const double sim_time)
{
    const auto start_tp = std::chrono::system_clock::from_time_t(start_time);
    const auto duration = std::chrono::duration_cast<std::chrono::system_clock::duration>(
        std::chrono::duration<double>(sim_time));
    const std::time_t cur_time = std::chrono::system_clock::to_time_t(start_tp + duration);
    std::ostringstream ss;
    ss << std::put_time(std::localtime(&cur_time), "%d-%b-%Y");
    return ss.str();
}

} // anonymous namespace

namespace Opm {

template<typename Scalar, typename IndexTraits>
template<class RatioFunc>
bool WellTest<Scalar, IndexTraits>::
checkMaxRatioLimitWell(const SingleWellState<Scalar, IndexTraits>& ws,
                       const Scalar max_ratio_limit,
                       const RatioFunc& ratioFunc,
                       Scalar& well_ratio_value) const
{
    const int np = well_.numPhases();

    std::vector<Scalar> well_rates(np, 0.0);
    for (int p = 0; p < np; ++p) {
        well_rates[p] = ws.surface_rates[p];
    }

    well_ratio_value = ratioFunc(well_rates, well_.phaseUsage());
    return (well_ratio_value > max_ratio_limit);
}

template<typename Scalar, typename IndexTraits>
template<class RatioFunc>
void WellTest<Scalar, IndexTraits>::
checkMaxRatioLimitCompletions(const SingleWellState<Scalar, IndexTraits>& ws,
                              const Scalar max_ratio_limit,
                              const Scalar well_ratio_value,
                              const RatioFunc& ratioFunc,
                              const std::string& ratio_name,
                              const UnitSystem::measure ratio_measure,
                              RatioLimitCheckReport& report) const
{
    int worst_offending_completion = RatioLimitCheckReport::INVALIDCOMPLETION;

    // the maximum water cut value of the completions
    // it is used to identify the most offending completion
    Scalar max_ratio_completion = 0;
    const int np = well_.numPhases();

    const auto& perf_data = ws.perf_data;
    const auto& perf_phase_rates = perf_data.phase_rates;
    // look for the worst_offending_completion
    for (const auto& completion : well_.getCompletions()) {
        std::vector<Scalar> completion_rates(np, 0.0);

        // looping through the connections associated with the completion
        const std::vector<int>& conns = completion.second;
        for (const int c : conns) {
            for (int p = 0; p < np; ++p) {
                const Scalar connection_rate = perf_phase_rates[c * np + p];
                completion_rates[p] += connection_rate;
            }
        } // end of for (const int c : conns)

        const Scalar ratio_completion = ratioFunc(completion_rates, well_.phaseUsage());

        if (ratio_completion > max_ratio_completion) {
            worst_offending_completion = completion.first;
            max_ratio_completion = ratio_completion;
        }
    } // end of for (const auto& completion : completions_)

    const Scalar local_max_ratio_completion = max_ratio_completion;
    max_ratio_completion = well_.parallelWellInfo().communication().max(max_ratio_completion);
    const Scalar violation_extent = max_ratio_completion / max_ratio_limit;

    if (violation_extent > report.violation_extent) {
        if (well_.parallelWellInfo().communication().size() > 1) {
            // well is distributed, communicate the worst-offending completion
            if (local_max_ratio_completion != max_ratio_completion)
                worst_offending_completion = std::numeric_limits<int>::min();
            worst_offending_completion = well_.parallelWellInfo().communication().max(worst_offending_completion);
        }

        report.worst_offending_completion = worst_offending_completion;
        report.violation_extent = violation_extent;
        // Report the well-level ratio (WECON); the completion ratio only
        // selects the worst offender via violation_extent above.
        report.ratio_name = ratio_name;
        report.ratio_measure = ratio_measure;
        report.ratio_value = well_ratio_value;
        report.ratio_limit = max_ratio_limit;
    }
}

template<typename Scalar, typename IndexTraits>
void WellTest<Scalar, IndexTraits>::
checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                 const SingleWellState<Scalar, IndexTraits>& ws,
                 RatioLimitCheckReport& report) const
{
    // function to calculate gor based on rates
    auto gor = [](const std::vector<Scalar>& rates,
                  const PhaseUsageInfo<IndexTraits>& pu)
    {
        const Scalar oil_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
        const Scalar gas_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)];
        if (gas_rate <= 0.)
            return Scalar{0};
        else if (oil_rate <= 0.)
            return Scalar{1e30}; // big value to mark it as violated
        else
            return (gas_rate / oil_rate);
    };

    const Scalar max_gor_limit = econ_production_limits.maxGasOilRatio();
    assert(max_gor_limit > 0.);

    Scalar well_gor = 0.0;
    const bool gor_limit_violated = this->checkMaxRatioLimitWell(ws, max_gor_limit, gor, well_gor);

    if (gor_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_gor_limit, well_gor, gor,
                                            "gas-oil ratio",
                                            UnitSystem::measure::gas_oil_ratio,
                                            report);
    }
}

template<typename Scalar, typename IndexTraits>
void WellTest<Scalar, IndexTraits>::
checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                 const SingleWellState<Scalar, IndexTraits>& ws,
                 RatioLimitCheckReport& report) const
{
    // function to calculate wgr based on rates
    auto wgr = [](const std::vector<Scalar>& rates,
                  const PhaseUsageInfo<IndexTraits>& pu)
    {
        const Scalar water_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
        const Scalar gas_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)];
        if (water_rate <= 0.)
            return Scalar{0};
        else if (gas_rate <= 0.)
            return Scalar{1e30}; // big value to mark it as violated
        else
            return (water_rate / gas_rate);
    };

    const Scalar max_wgr_limit = econ_production_limits.maxWaterGasRatio();
    assert(max_wgr_limit > 0.);

    Scalar well_wgr = 0.0;
    const bool wgr_limit_violated = this->checkMaxRatioLimitWell(ws, max_wgr_limit, wgr, well_wgr);

    if (wgr_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_wgr_limit, well_wgr, wgr,
                                            "water-gas ratio",
                                            UnitSystem::measure::water_gas_ratio,
                                            report);
    }
}

template<typename Scalar, typename IndexTraits>
void WellTest<Scalar, IndexTraits>::
checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                      const SingleWellState<Scalar, IndexTraits>& ws,
                      RatioLimitCheckReport& report) const
{
    // function to calculate water cut based on rates
    auto waterCut = [](const std::vector<Scalar>& rates,
                       const PhaseUsageInfo<IndexTraits>& pu)
    {
        const Scalar oil_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
        const Scalar water_rate = -rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
        const Scalar liquid_rate = oil_rate + water_rate;
        if (liquid_rate <= 0.)
            return Scalar{0};
        else if (water_rate < 0)
            return Scalar{0};
        else if (oil_rate < 0)
            return Scalar{1};
        else
            return (water_rate / liquid_rate);

    };

    const Scalar max_water_cut_limit = econ_production_limits.maxWaterCut();
    assert(max_water_cut_limit > 0.);

    Scalar well_water_cut = 0.0;
    const bool watercut_limit_violated =
        this->checkMaxRatioLimitWell(ws, max_water_cut_limit, waterCut, well_water_cut);

    if (watercut_limit_violated) {
        report.ratio_limit_violated = true;
        this->checkMaxRatioLimitCompletions(ws, max_water_cut_limit, well_water_cut, waterCut,
                                            "water cut",
                                            UnitSystem::measure::water_cut,
                                            report);
    }
}

template<typename Scalar, typename IndexTraits>
bool WellTest<Scalar, IndexTraits>::
checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                    const std::vector<Scalar>& rates_or_potentials,
                    DeferredLogger& deferred_logger) const
{
    const auto& pu = well_.phaseUsage();
    if (econ_production_limits.onMinOilRate()) {
        const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        const Scalar oil_rate = rates_or_potentials[oil_pos];
        const Scalar min_oil_rate = econ_production_limits.minOilRate();
        if (std::abs(oil_rate) < min_oil_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinGasRate() ) {
        const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        const Scalar gas_rate = rates_or_potentials[gas_pos];
        const Scalar min_gas_rate = econ_production_limits.minGasRate();
        if (std::abs(gas_rate) < min_gas_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinLiquidRate() ) {
        const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        const Scalar oil_rate = rates_or_potentials[oil_pos];
        const Scalar water_rate = rates_or_potentials[water_pos];
        const Scalar liquid_rate = oil_rate + water_rate;
        const Scalar min_liquid_rate = econ_production_limits.minLiquidRate();
        if (std::abs(liquid_rate) < min_liquid_rate) {
            return true;
        }
    }

    if (econ_production_limits.onMinReservoirFluidRate()) {
        deferred_logger.warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
    }

    return false;
}

template<typename Scalar, typename IndexTraits>
typename WellTest<Scalar, IndexTraits>::RatioLimitCheckReport
WellTest<Scalar, IndexTraits>::
checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                     const SingleWellState<Scalar, IndexTraits>& ws,
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

template<typename Scalar, typename IndexTraits>
void WellTest<Scalar, IndexTraits>::
updateWellTestStateEconomic(const SingleWellState<Scalar, IndexTraits>& ws,
                            const double simulation_time,
                            const bool write_message_to_opmlog,
                            WellTestState& well_test_state,
                            const bool zero_group_target,
                            const UnitSystem& unit_system,
                            const std::time_t start_time,
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
            if (!zero_group_target) {
                rate_limit_violated
                    = this->checkRateEconLimits(econ_production_limits, ws.surface_rates, deferred_logger);
            }
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

        // Shared "at time ..." and ratio-violation clauses for the messages
        // below; only built when a message will actually be logged.
        std::string when;
        std::string reason;
        if (write_message_to_opmlog) {
            when = fmt::format(
                "at time {:.2f} {} (date = {})",
                unit_system.from_si(UnitSystem::measure::time, simulation_time),
                unit_system.name(UnitSystem::measure::time),
                dateString(start_time, simulation_time));

            const std::string ratio_unit = unit_system.name(ratio_report.ratio_measure);
            const std::string unit_suffix = ratio_unit.empty() ? std::string{}
                                                               : " " + ratio_unit;
            reason = fmt::format(
                "{} {:.4e}{} exceeds the limit {:.4e}{}",
                ratio_report.ratio_name,
                unit_system.from_si(ratio_report.ratio_measure, ratio_report.ratio_value), unit_suffix,
                unit_system.from_si(ratio_report.ratio_measure, ratio_report.ratio_limit), unit_suffix);
        }

        switch (workover) {
        case WellEconProductionLimits::EconWorkover::CON:
        case WellEconProductionLimits::EconWorkover::CONP:
            {
                // CON  : shut the worst-offending connection/completion only.
                // +CON : shut the worst-offending connection/completion and all
                //        connections below it in the wellbore.
                const bool close_connections_below =
                    (workover == WellEconProductionLimits::EconWorkover::CONP);
                this->closeOffendingCompletion(ratio_report.worst_offending_completion,
                                               close_connections_below,
                                               simulation_time,
                                               write_message_to_opmlog,
                                               well_test_state,
                                               when, reason,
                                               deferred_logger);
                break;
            }
        case WellEconProductionLimits::EconWorkover::WELL:
            {
            well_test_state.close_well(well_.name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
            if (write_message_to_opmlog) {
                const std::string action = well_.wellEcl().getAutomaticShutIn() ? "shut" : "stopped";
                const std::string sep = this->message_separator();
                deferred_logger.info(
                    fmt::format("{}\nWell {} will be {} {},\nBecause {}.\n{}",
                                sep, well_.name(), action, when, reason, sep));
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

template<typename Scalar, typename IndexTraits>
void WellTest<Scalar, IndexTraits>::
updateWellTestStateCECON(const SingleWellState<Scalar, IndexTraits>& ws,
                         const double simulation_time,
                         const bool write_message_to_opmlog,
                         WellTestState& well_test_state,
                         const UnitSystem& unit_system,
                         const std::time_t start_time,
                         DeferredLogger& deferred_logger) const
{
    if (well_.isInjector())
        return;

    // The well may already have been closed by the WECON checks in
    // updateWellTestStateEconomic().
    if (well_test_state.well_is_closed(well_.name()))
        return;

    const auto& completions = well_.getCompletions(); // complnum -> [rank-local perf indices]
    const int np = well_.numPhases();
    const auto& perf_phase_rates = ws.perf_data.phase_rates;
    const auto& pu = well_.phaseUsage();
    const auto& comm = well_.parallelWellInfo().communication();

    std::vector<Scalar> conn_rates(np);
    for (const auto& connection : well_.wellEcl().getConnections()) {
        if (!connection.hasEconLimits())
            continue;

        // The skip conditions below are based on deck data and the well test
        // state, which are identical on all ranks sharing a distributed well,
        // so every rank participates in the rate summation for the same
        // sequence of completions.
        if (connection.state() != Connection::State::OPEN)
            continue;

        const int complnum = connection.complnum();
        if (well_test_state.completion_is_closed(well_.name(), complnum))
            continue;

        // Note: minimum oil/gas rate limits (CECON items 11 and 12) are not
        // supported; their use is flagged at deck validation time.
        const auto& limits = connection.econLimits();

        // Sum rates over all perforations belonging to this completion. For a
        // distributed well each rank contributes its local perforations and the
        // result is summed over the well communicator, so that all ranks reach
        // the same closing decisions below.
        std::fill(conn_rates.begin(), conn_rates.end(), Scalar{0});
        if (const auto it = completions.find(complnum); it != completions.end()) {
            for (const int perf : it->second) {
                for (int p = 0; p < np; ++p) {
                    conn_rates[p] += perf_phase_rates[perf * np + p];
                }
            }
        }
        if (comm.size() > 1) {
            comm.sum(conn_rates.data(), np);
        }

        // Production rates are stored as negative values; negate for the checks below.
        const Scalar oil_rate = pu.phaseIsActive(IndexTraits::oilPhaseIdx)
            ? -conn_rates[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)]
            : Scalar{0};
        const Scalar gas_rate = pu.phaseIsActive(IndexTraits::gasPhaseIdx)
            ? -conn_rates[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)]
            : Scalar{0};
        const Scalar water_rate = pu.phaseIsActive(IndexTraits::waterPhaseIdx)
            ? -conn_rates[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)]
            : Scalar{0};

        // Record which ratio limit was violated together with its actual value
        // and limit, so the closing message can report the offending quantity.
        bool violated = false;
        const char* ratio_name = nullptr;
        UnitSystem::measure ratio_measure = UnitSystem::measure::identity;
        Scalar ratio_value = Scalar{0};
        Scalar ratio_limit = Scalar{0};

        if (!violated && limits.onMaxWaterCut()) {
            const Scalar liquid = oil_rate + water_rate;
            if (liquid > Scalar{0}) {
                const Scalar wcut = (oil_rate >= Scalar{0})
                    ? water_rate / liquid
                    : Scalar{1};
                if (wcut > limits.max_water_cut) {
                    violated = true;
                    ratio_name = "water cut";
                    ratio_measure = UnitSystem::measure::water_cut;
                    ratio_value = wcut;
                    ratio_limit = limits.max_water_cut;
                }
            }
        }

        if (!violated && limits.onMaxGasOilRatio()) {
            if (gas_rate > Scalar{0}) {
                const bool no_oil = (oil_rate <= Scalar{0});
                if (no_oil || (gas_rate / oil_rate > limits.max_gas_oil_ratio)) {
                    violated = true;
                    ratio_name = "gas-oil ratio";
                    ratio_measure = UnitSystem::measure::gas_oil_ratio;
                    ratio_value = no_oil ? std::numeric_limits<Scalar>::infinity()
                                         : gas_rate / oil_rate;
                    ratio_limit = limits.max_gas_oil_ratio;
                }
            }
        }

        if (!violated && limits.onMaxWaterGasRatio()) {
            if (water_rate > Scalar{0}) {
                const bool no_gas = (gas_rate <= Scalar{0});
                if (no_gas || (water_rate / gas_rate > limits.max_water_gas_ratio)) {
                    violated = true;
                    ratio_name = "water-gas ratio";
                    ratio_measure = UnitSystem::measure::water_gas_ratio;
                    ratio_value = no_gas ? std::numeric_limits<Scalar>::infinity()
                                         : water_rate / gas_rate;
                    ratio_limit = limits.max_water_gas_ratio;
                }
            }
        }

        if (!violated)
            continue;

        // Build the "at time ... (date = ...)" and ratio-violation clauses that
        // are shared by all CECON workover messages below.
        const std::string when = fmt::format(
            "at time {:.2f} {} (date = {})",
            unit_system.from_si(UnitSystem::measure::time, simulation_time),
            unit_system.name(UnitSystem::measure::time),
            ceconDateString(start_time, simulation_time));

        const std::string ratio_unit = unit_system.name(ratio_measure);
        const std::string unit_suffix = ratio_unit.empty() ? std::string{}
                                                           : " " + ratio_unit;
        const std::string reason = fmt::format(
            "{} {:.4e}{} exceeds the limit {:.4e}{}",
            ratio_name,
            unit_system.from_si(ratio_measure, ratio_value), unit_suffix,
            unit_system.from_si(ratio_measure, ratio_limit), unit_suffix);

        switch (limits.workover) {
        case ConnectionEconLimits::EconWorkover::CON:
        case ConnectionEconLimits::EconWorkover::CONP:
            {
                const bool close_below =
                    (limits.workover == ConnectionEconLimits::EconWorkover::CONP);
                this->closeOffendingCompletion(complnum, close_below,
                                                simulation_time, write_message_to_opmlog,
                                                well_test_state, when, reason,
                                                deferred_logger);
                if (well_test_state.well_is_closed(well_.name()))
                    return;
                break;
            }
        case ConnectionEconLimits::EconWorkover::WELL:
            well_test_state.close_well(well_.name(), WellTestConfig::Reason::ECONOMIC,
                                       simulation_time);
            if (write_message_to_opmlog) {
                const std::string action = well_.wellEcl().getAutomaticShutIn() ? "shut" : "stopped";
                const std::string sep = this->message_separator();
                deferred_logger.info(
                    fmt::format("{}\nWell {} will be {} {},\nBecause {} {}.\n{}",
                                sep, well_.name(), action, when,
                                this->completionDescriptor(complnum), reason, sep));
            }
            return;
        case ConnectionEconLimits::EconWorkover::PLUG:
            deferred_logger.warning("NOT_SUPPORTING_CECON_PLUG",
                                    "CECON PLUG workover is not supported yet");
            break;
        }
    }
}

template<typename Scalar, typename IndexTraits>
void WellTest<Scalar, IndexTraits>::
closeOffendingCompletion(const int offending_completion,
                         const bool close_connections_below,
                         const double simulation_time,
                         const bool write_message_to_opmlog,
                         WellTestState& well_test_state,
                         const std::string& when,
                         const std::string& reason,
                         DeferredLogger& deferred_logger) const
{
    // complnum is always >= 1; a non-positive value would be a bug.
    assert(offending_completion > 0);

    // Sentinel for "no offending completion"; normally filtered by the caller.
    if (offending_completion == RatioLimitCheckReport::INVALIDCOMPLETION) {
        return;
    }

    // The wellbore ordering of the connections is the one given by COMPORD and
    // reflected in the order of the connections returned by getConnections().
    const auto& connections = well_.wellEcl().getConnections();

    // Build the set of completions to close. For CON it is just the
    // worst-offending completion; for +CON it additionally contains every
    // completion located below the worst-offending one in the wellbore.
    // "Below" means further from the wellhead according to the connection
    // ordering. A set is used because a completion spanning several connections
    // would otherwise be inserted multiple times.
    std::unordered_set<int> completions_to_close{offending_completion};
    if (close_connections_below) {
        bool below_worst_offender = false;
        for (const auto& connection : connections) {
            if (!below_worst_offender &&
                connection.complnum() == offending_completion) {
                below_worst_offender = true;
                continue;
            }
            if (below_worst_offender) {
                completions_to_close.insert(connection.complnum());
            }
        }
    }

    for (const int completion : completions_to_close) {
        well_test_state.close_completion(well_.name(), completion, simulation_time);
    }

    if (write_message_to_opmlog) {
        // "and all below" is only appended for the +CON workover, which also
        // closes every connection below the offending completion.
        const std::string subject = close_connections_below
            ? fmt::format("{}, and all below, in Well {}",
                          this->completionDescriptor(offending_completion), well_.name())
            : fmt::format("{} in Well {}",
                          this->completionDescriptor(offending_completion), well_.name());
        const std::string sep = this->message_separator();
        // WECON-driven: report "the well" ratio (CECON would say "its ...", will be addressed by CECON development).
        deferred_logger.info(
            fmt::format("{}\n{} will be closed {},\nBecause the well {}.\n{}",
                        sep, subject, when, reason, sep));
    }

    bool allCompletionsClosed = true;
    for (const auto& connection : connections) {
        if (connection.state() == Connection::State::OPEN
            && !well_test_state.completion_is_closed(well_.name(), connection.complnum())) {
            allCompletionsClosed = false;
        }
    }

    if (allCompletionsClosed) {
        well_test_state.close_well(well_.name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
        // if all the completion/connections are closed, the well can only be SHUT
        if (write_message_to_opmlog) {
            const std::string sep = this->message_separator();
            deferred_logger.info(
                fmt::format("{}\nWell {} will be shut {},\n"
                            "Because all its completions are closed.\n{}",
                            sep, well_.name(), when, sep));
        }
    }
}

template<typename Scalar, typename IndexTraits>
std::string WellTest<Scalar, IndexTraits>::
completionDescriptor(const int complnum) const
{
    const auto& connections = well_.wellEcl().getConnections();
    const Connection* first = nullptr;
    int num_connections = 0;
    for (const auto& connection : connections) {
        if (connection.complnum() == complnum) {
            if (first == nullptr) {
                first = &connection;
            }
            ++num_connections;
        }
    }

    // A completion that spans several connections/blocks is identified by its
    // number only; a single-connection completion also reports its block.
    if (num_connections == 1 && first != nullptr) {
        return fmt::format("Completion {} - block ({}, {}, {})", complnum,
                           first->getI() + 1,
                           first->getJ() + 1,
                           first->getK() + 1);
    }
    return fmt::format("Completion {}", complnum);
}

template<typename Scalar, typename IndexTraits>
void WellTest<Scalar, IndexTraits>::
updateWellTestStatePhysical(const double simulation_time,
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

template class WellTest<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class WellTest<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
