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
{
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

template class StandardWellGeneric<double>;

}
