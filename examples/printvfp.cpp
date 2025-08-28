/*
  Copyright 2020 Equinor.

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

#include <opm/common/utility/TimeService.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Opm;

namespace {

struct Setup
{
    using IndexTraits = BlackOilDefaultFluidSystemIndices;
    using PhaseUsage = PhaseUsageInfo<IndexTraits>;
    explicit Setup(const std::string& file)
        : Setup { Parser{}.parseFile(file) }
    {}

    explicit Setup(const Deck& deck)
        : ecl_state  { std::make_unique<EclipseState>(deck) }
        , schedule   { std::make_unique<Schedule>(deck, *ecl_state, std::make_shared<Python>()) }
    {
        const int step = 0;
        const auto& sched_state = (*this->schedule)[step];

        PhaseUsage pu;
        pu.initFromPhases(ecl_state->runspec().phases());

        this->well_state = std::make_unique<WellState<double, IndexTraits>>(pu);

        this->vfp_properties = std::make_unique<VFPProperties<double, IndexTraits>>
            (sched_state.vfpinj(), sched_state.vfpprod(), *well_state);
    }

    std::unique_ptr<EclipseState> ecl_state;
    std::unique_ptr<Schedule> schedule;
    std::unique_ptr<WellState<double, IndexTraits>> well_state;
    std::unique_ptr<VFPProperties<double, IndexTraits>> vfp_properties;
};


double computeBhp(const VFPProdTable& table,
                  const double flo,
                  const double thp,
                  const double wfr,
                  const double gfr,
                  const double alq)
{

    // First, find the values to interpolate between.
    // Assuming positive flo here!
    assert(flo > 0.0);

    const auto flo_i = VFPHelpers<double>::findInterpData(flo, table.getFloAxis());
    const auto thp_i = VFPHelpers<double>::findInterpData(thp, table.getTHPAxis()); // assume constant
    const auto wfr_i = VFPHelpers<double>::findInterpData(wfr, table.getWFRAxis());
    const auto gfr_i = VFPHelpers<double>::findInterpData(gfr, table.getGFRAxis());
    const auto alq_i = VFPHelpers<double>::findInterpData(alq, table.getALQAxis()); // assume constant

    return VFPHelpers<double>::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
}

} // Anonymous namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        return EXIT_FAILURE;
    }

    const Setup setup(argv[1]);

//     const int table_id = 1;
    const int table_id = 4;
    const double wct = 0.0;
//    const double gor = 35.2743;
    const double gor = 0.0;
    const double alq = 0.0;
    const int n = 51;
    const double m3pd = unit::cubic(unit::meter)/unit::day;
    const double rate_min = 20.0 * m3pd;
    const double rate_max = 2000.0 * m3pd;
    // const double thp = 32.1744 * unit::barsa;
    // const double thp = 10.0 * unit::barsa;
    const double thp_min = 10.0 * unit::barsa;
    const double thp_max = 35.0 * unit::barsa;
    std::vector<double> rates(n);
    std::vector<double> thps(n);
    for (int ii = 0; ii < n; ++ii) {
        const double q = double(ii) / double(n-1);
        rates[ii] = (1.0 - q) * rate_min + q * rate_max;
        thps[ii] = (1.0 - q) * thp_min + q * thp_max;
    }

    const VFPProdTable& table = setup.vfp_properties->getProd()->getTable(table_id);
    std::cout.precision(12);
    for (double rate : rates) {
        for (double thp : thps) {
            const double bhp = computeBhp(table, rate, thp, wct, gor, alq);
            std::cout //<< std::setw(18) << unit::convert::to(rate, m3pd)
                      //<< std::setw(18) << unit::convert::to(thp, unit::barsa)
                      << std::setw(18) << unit::convert::to(bhp, unit::barsa) << '\n';
        }
    }

    return EXIT_SUCCESS;
}
