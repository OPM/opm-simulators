/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

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
#include <opm/simulators/flow/NonlinearSolverEbos.hpp>

#include <cmath>

namespace Opm {
namespace detail {

void detectOscillations(const std::vector<std::vector<double>>& residualHistory,
                        const int it, const int numPhases, const double relaxRelTol,
                        bool& oscillate, bool& stagnate)
{
    // The detection of oscillation in two primary variable results in the report of the detection
    // of oscillation for the solver.
    // Only the saturations are used for oscillation detection for the black oil model.
    // Stagnate is not used for any treatment here.

    if (it < 2) {
        oscillate = false;
        stagnate = false;
        return;
    }

    stagnate = true;
    int oscillatePhase = 0;
    const auto& F0 = residualHistory[it];
    const auto& F1 = residualHistory[it - 1];
    const auto& F2 = residualHistory[it - 2];
    for (int p = 0; p < numPhases; ++p) {
        const double d1 = std::abs((F0[p] - F2[p]) / F0[p]);
        const double d2 = std::abs((F0[p] - F1[p]) / F0[p]);

        oscillatePhase += (d1 < relaxRelTol) && (relaxRelTol < d2);

        // Process is 'stagnate' unless at least one phase
        // exhibits significant residual change.
        stagnate = (stagnate && !(std::abs((F1[p] - F2[p]) / F2[p]) > 1.0e-3));
    }

    oscillate = (oscillatePhase > 1);
}

} // namespace detail
} // namespace Opm
