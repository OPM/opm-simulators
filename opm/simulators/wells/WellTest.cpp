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

#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

namespace Opm
{

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

} // namespace Opm
