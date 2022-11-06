/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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
#include <opm/simulators/wells/StandardWellGeneric.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>
#include <stdexcept>

namespace Opm
{

template<class Scalar>
StandardWellGeneric<Scalar>::
StandardWellGeneric(const WellInterfaceGeneric& baseif)
    : baseif_(baseif)
    , perf_densities_(baseif_.numPerfs())
    , perf_pressure_diffs_(baseif_.numPerfs())
    , parallelB_(duneB_, baseif_.parallelWellInfo())
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise);
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}


template<class Scalar>
double
StandardWellGeneric<Scalar>::
relaxationFactorRate(const std::vector<double>& primary_variables,
                     const BVectorWell& dwells)
{
    double relaxation_factor = 1.0;
    static constexpr int WQTotal = 0;

    // For injector, we only check the total rates to avoid sign change of rates
    const double original_total_rate = primary_variables[WQTotal];
    const double newton_update = dwells[0][WQTotal];
    const double possible_update_total_rate = primary_variables[WQTotal] - newton_update;

    // 0.8 here is a experimental value, which remains to be optimized
    // if the original rate is zero or possible_update_total_rate is zero, relaxation_factor will
    // always be 1.0, more thoughts might be needed.
    if (original_total_rate * possible_update_total_rate < 0.) { // sign changed
        relaxation_factor = std::abs(original_total_rate / newton_update) * 0.8;
    }

    assert(relaxation_factor >= 0.0 && relaxation_factor <= 1.0);

    return relaxation_factor;
}

template<class Scalar>
double
StandardWellGeneric<Scalar>::
relaxationFactorFraction(const double old_value,
                         const double dx)
{
    assert(old_value >= 0. && old_value <= 1.0);

    double relaxation_factor = 1.;

    // updated values without relaxation factor
    const double possible_updated_value = old_value - dx;

    // 0.95 is an experimental value remains to be optimized
    if (possible_updated_value < 0.0) {
        relaxation_factor = std::abs(old_value / dx) * 0.95;
    } else if (possible_updated_value > 1.0) {
        relaxation_factor = std::abs((1. - old_value) / dx) * 0.95;
    }
    // if possible_updated_value is between 0. and 1.0, then relaxation_factor
    // remains to be one

    assert(relaxation_factor >= 0. && relaxation_factor <= 1.);

    return relaxation_factor;
}

template<class Scalar>
void
StandardWellGeneric<Scalar>::
computeConnectionPressureDelta()
{
    // Algorithm:

    // We'll assume the perforations are given in order from top to
    // bottom for each well.  By top and bottom we do not necessarily
    // mean in a geometric sense (depth), but in a topological sense:
    // the 'top' perforation is nearest to the surface topologically.
    // Our goal is to compute a pressure delta for each perforation.

    // 1. Compute pressure differences between perforations.
    //    dp_perf will contain the pressure difference between a
    //    perforation and the one above it, except for the first
    //    perforation for each well, for which it will be the
    //    difference to the reference (bhp) depth.

    const int nperf = baseif_.numPerfs();
    perf_pressure_diffs_.resize(nperf, 0.0);
    auto z_above = baseif_.parallelWellInfo().communicateAboveValues(baseif_.refDepth(), baseif_.perfDepth());

    for (int perf = 0; perf < nperf; ++perf) {
        const double dz = baseif_.perfDepth()[perf] - z_above[perf];
        perf_pressure_diffs_[perf] = dz * perf_densities_[perf] * baseif_.gravity();
    }

    // 2. Compute pressure differences to the reference point (bhp) by
    //    accumulating the already computed adjacent pressure
    //    differences, storing the result in dp_perf.
    //    This accumulation must be done per well.
    const auto beg = perf_pressure_diffs_.begin();
    const auto end = perf_pressure_diffs_.end();
    baseif_.parallelWellInfo().partialSumPerfValues(beg, end);
}



template<class Scalar>
unsigned int
StandardWellGeneric<Scalar>::
getNumBlocks() const
{
    return duneB_.nonzeroes();
}

template class StandardWellGeneric<double>;

}
