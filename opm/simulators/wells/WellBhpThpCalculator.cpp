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
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

namespace Opm
{

bool
WellBhpThpCalculator::wellHasTHPConstraints(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto controls = well_ecl.injectionControls(summaryState);
        if (controls.hasControl(Well::InjectorCMode::THP))
            return true;
    }

    if (well_ecl.isProducer()) {
        const auto controls = well_ecl.productionControls(summaryState);
        if (controls.hasControl(Well::ProducerCMode::THP))
            return true;
    }

    return false;
}

double WellBhpThpCalculator::getTHPConstraint(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto& controls = well_ecl.injectionControls(summaryState);
        return controls.thp_limit;
    }

    if (well_ecl.isProducer( )) {
        const auto& controls = well_ecl.productionControls(summaryState);
        return controls.thp_limit;
    }

    return 0.0;
}

double WellBhpThpCalculator::mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto& controls = well_ecl.injectionControls(summaryState);
        return controls.bhp_limit;
    }

    if (well_ecl.isProducer( )) {
        const auto& controls = well_ecl.productionControls(summaryState);
        return controls.bhp_limit;
    }

    return 0.0;
}

} // namespace Opm
